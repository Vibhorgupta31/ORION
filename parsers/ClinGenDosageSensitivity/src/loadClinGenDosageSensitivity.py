import os
import enum
import gzip
from Common.extractor import Extractor
from Common.loader_interface import SourceDataLoader
from biolink_constants import PRIMARY_KNOWLEDGE_SOURCE, NODE_TYPES, SEQUENCE_VARIANT
from Common.utils import GetData
from datetime import date


# For parsing tsv files
# HI: HaploInsufficiency ; TS : TriploSensitivity
class ClinGenDosageSensitivityCOLS(enum.IntEnum):
    REGION = 0
    GENE = 1
    HI_DISEASE = -2
    TS_DISEASE = -1
    HI_SCORE = 4
    HI_DESCRIPTION = 5
    TS_SCORE = 12
    TS_DESCRIPTION = 13


##############
# Class: ClinGen Dosage Sensitivity  data loader
# Desc: Class that loads/parses the ClinGen Dosage Sensitivity  data.
##############
class ClinGenDosageSensitivityLoader(SourceDataLoader):
    source_id: str = "ClinGenDosageSensitivity"
    provenance_id: str = (
        "infores:clingen"  # Need to figure out this, can only be filled from one of the values from https://github.com/biolink/biolink-model/blob/master/infores_catalog.yaml
    )
    # increment parsing_version whenever changes are made to the parser that would result in changes to parsing output
    parsing_version: str = "v1.0"

    def __init__(self, test_mode: bool = False, source_data_dir: str = None):
        """
        :param test_mode - sets the run into test mode
        :param source_data_dir - the specific storage directory to save files in
        """
        super().__init__(test_mode=test_mode, source_data_dir=source_data_dir)
        self.cligen_dosage_sensitivity_url = "http://ftp.clinicalgenome.org/"
        self.cligen_dosage_sensitivity_gene_file = (
            "ClinGen_gene_curation_list_GRCh38.tsv"
        )
        self.clingen_dosage_sensitivity_region_file = (
            "ClinGen_region_curation_list_GRCh38.tsv"
        )
        self.data_files = [
            self.cligen_dosage_sensitivity_gene_file,
            self.clingen_dosage_sensitivity_region_file,
        ]

    def get_latest_source_version(self) -> str:
        # if possible go to the source and retrieve a string that is the latest version of the source data
        latest_version = date.today().strftime("%Y%m%d")
        return latest_version

    def get_data(self) -> bool:
        # get_data is responsible for fetching the files in self.data_files and saving them to self.data_path
        data_puller = GetData()
        for source in self.data_files:
            source_data_url = f"{self.cligen_dosage_sensitivity_url}{source}"
            data_puller.pull_via_http(source_data_url, self.data_path)
        return True

    def process_dosage_sensitivity_file(self, extractor, filename: str, subject_extractor, predicate: str):
        with open(filename, "rt") as fp:
            extractor.csv_extract(
                fp,
                subject_extractor,
                # subject id
                lambda line: (
                    "MONDO:0700096"
                    if ((line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value] == "") and (
                            line[ClinGenDosageSensitivityCOLS.TS_DISEASE.value] == ""))
                    else (line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value] if (
                            line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value] != "") else line[
                        ClinGenDosageSensitivityCOLS.TS_DISEASE.value])),
                # object id
                lambda line: predicate,  # predicate extractor
                lambda line: {},  # subject properties
                lambda line: {},  # object properties
                lambda line: dict(
                    {
                        PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                        "Haploinsufficiency Description": line[
                            ClinGenDosageSensitivityCOLS.HI_DESCRIPTION.value
                        ],
                        "Triplosensitivity Description": line[
                            ClinGenDosageSensitivityCOLS.TS_DESCRIPTION.value
                        ],
                    },
                    **get_edge_properties(
                        line[ClinGenDosageSensitivityCOLS.HI_SCORE.value],
                        line[ClinGenDosageSensitivityCOLS.TS_SCORE.value],
                        line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value],
                        line[ClinGenDosageSensitivityCOLS.TS_DISEASE.value],
                    ),
                ),
                # edge properties
                comment_character="#",
                delim="\t",
                has_header_row=True,
            )

    def parse_data(self) -> dict:
        """
        Parses the data file for graph nodes/edges

        :return: ret_val: load_metadata
        """
        # This is a made up example of how one might extract nodes and edges from a tsv file
        # In this case it's taking the subject ID from column 1 and the object ID from column 3,
        # prepending them with a curie prefix. The predicate comes from column 2. The value in column 4
        # is set as a property on the edge.
        extractor = Extractor(file_writer=self.output_file_writer)
        dosage_sensitivity_gene_file: str = os.path.join(
            self.data_path, self.cligen_dosage_sensitivity_gene_file
        )
        dosage_sensitivity_region_file: str = os.path.join(
            self.data_path, self.clingen_dosage_sensitivity_region_file
        )
        
        self.process_dosage_sensitivity_file(extractor, dosage_sensitivity_gene_file,
                                             lambda line: f"NCBIGene:{line[ClinGenDosageSensitivityCOLS.GENE.value]}",
                                             "gene associated with condition")
        self.process_dosage_sensitivity_file(extractor, dosage_sensitivity_region_file,
                                            lambda line: f"ISCARegion:{line[ClinGenDosageSensitivityCOLS.REGION.value]}",
                                            "region associated with condition")

        return extractor.load_metadata


# Created a common function to take the score value and return attributes,
# this may have problems with TS and HI score where the mode of inheritance is captured, followed ClinGen criteria for converting scores
# Scoring values are from https://clinicalgenome.org/site/assets/files/6428/dosage_sop-scoring-1.pdf
def get_edge_properties(hi_score, ts_score, hi_mondo_id, ts_mondo_id):
    if hi_mondo_id != "" or ts_mondo_id != "":
        try:
            hi_score = int(hi_score)
            ts_score = int(ts_score)
        except ValueError:
            return {"Status": "Not yet evaluated"}
        if (hi_score in (0, 1, 2, 3)) or (ts_score in (0, 1, 2, 3)):
            return {"negated": False}
        elif hi_score == 30 or ts_score == 30:
            return {"negated": False, "mode_of_inheritence": "Autosomal Recessive"}
        elif hi_score == 40 or ts_score == 40:  # another condition for negation
            return {"negated": True}
    else:  # Negate only if no HI disease and TS disease is available for the gene
        return {"negated": True}