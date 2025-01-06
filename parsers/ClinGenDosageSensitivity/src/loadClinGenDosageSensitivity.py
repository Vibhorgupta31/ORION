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
class ClinGenDosageSensitivityGeneCOLS(enum.IntEnum):
    GENE = 1
    HI_DISEASE_GENE = -2
    TS_DISEASE_GENE = -1
    HI_SCORE_GENE = 4
    HI_DESCRIPTION_GENE = 5
    TS_SCORE_GENE = 12
    TS_DESCRIPTION_GENE = 13


class ClinGenDosageSensitivityRegionCOLS(enum.IntEnum):
    REGION = 0
    HI_DISEASE_REGION = -2
    TS_DISEASE_REGION = -1
    HI_SCORE_REGION = 4
    HI_DESCRIPTION_REGION = 5
    TS_SCORE_REGION = 12
    TS_DESCRIPTION_REGION = 13


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
        with open(dosage_sensitivity_gene_file, "rt") as fp:
            extractor.csv_extract(
                fp,
                lambda line: f"NCBIGene:{line[ClinGenDosageSensitivityGeneCOLS.GENE.value]}",
                # subject id
                lambda line: (
                    "MONDO:0700096"
                    if ((line[ClinGenDosageSensitivityGeneCOLS.HI_DISEASE_GENE.value] == "") and (line[ClinGenDosageSensitivityGeneCOLS.TS_DISEASE_GENE.value]== ""))
                    else (line[ClinGenDosageSensitivityGeneCOLS.HI_DISEASE_GENE.value] if (line[ClinGenDosageSensitivityGeneCOLS.HI_DISEASE_GENE.value]!= "")else line[ClinGenDosageSensitivityGeneCOLS.TS_DISEASE_GENE.value])),
                # object id
                lambda line: "gene associated with condition",  # predicate extractor
                lambda line: {},  # subject properties
                lambda line: {},  # object properties
                lambda line: dict(
                    {
                        PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                        "Haploinsufficiency Description": line[
                            ClinGenDosageSensitivityGeneCOLS.HI_DESCRIPTION_GENE.value
                        ],
                        "Triplosensitivity Description": line[
                            ClinGenDosageSensitivityGeneCOLS.TS_DESCRIPTION_GENE.value
                        ],
                    },
                    **get_edge_properties(
                        line[ClinGenDosageSensitivityGeneCOLS.HI_SCORE_GENE.value],
                        line[ClinGenDosageSensitivityGeneCOLS.TS_SCORE_GENE.value],
                        line[ClinGenDosageSensitivityGeneCOLS.HI_DISEASE_GENE.value],
                        line[ClinGenDosageSensitivityGeneCOLS.TS_DISEASE_GENE.value],
                    ),
                ),
                # edge properties
                comment_character="#",
                delim="\t",
                has_header_row=True,
            )

        with open(dosage_sensitivity_region_file, "rt") as fp:
            extractor.csv_extract(
                fp,
                lambda line: f"ISCARegion:{line[ClinGenDosageSensitivityRegionCOLS.REGION.value]}",
                # subject id
                lambda line: (
                    "MONDO:0700096" if ((line[ClinGenDosageSensitivityRegionCOLS.HI_DISEASE_REGION.value] == "" ) and ( line[ClinGenDosageSensitivityRegionCOLS.TS_DISEASE_REGION.value]== ""))
                    else (line[ClinGenDosageSensitivityRegionCOLS.HI_DISEASE_REGION.value] if (line[ClinGenDosageSensitivityRegionCOLS.HI_DISEASE_REGION.value]!= "")else line[ClinGenDosageSensitivityRegionCOLS.TS_DISEASE_REGION.value])),
                # object id
                lambda line: "region associated with condition",  # predicate extractor
                lambda line: {},  # subject properties
                lambda line: {},  # object properties
                lambda line: dict(
                    {
                        PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                        "Haploinsufficiency Description": line[
                            ClinGenDosageSensitivityRegionCOLS.HI_DESCRIPTION_REGION.value
                        ],
                        "Triplosensitivity Description": line[
                            ClinGenDosageSensitivityRegionCOLS.TS_DESCRIPTION_REGION.value
                        ],
                    },
                    **get_edge_properties(
                        line[ClinGenDosageSensitivityRegionCOLS.HI_SCORE_REGION.value],
                        line[ClinGenDosageSensitivityRegionCOLS.TS_SCORE_REGION.value],
                        line[
                            ClinGenDosageSensitivityRegionCOLS.HI_DISEASE_REGION.value
                        ],
                        line[
                            ClinGenDosageSensitivityRegionCOLS.TS_DISEASE_REGION.value
                        ],
                    ),
                ),
                # edge properties
                comment_character="#",
                delim="\t",
                has_header_row=True,
            )
        return extractor.load_metadata


# Created a common function to take the score value and return attributes,
# this may have problems with TS and HI score where the mode of inheritance is captured, followed ClinGen criteria for converting scores
def get_edge_properties(hi_score, ts_score, hi_mondo_id, ts_mondo_id):
    if hi_mondo_id != "" or ts_mondo_id != "":
        try:
            hi_score = int(hi_score)
            ts_score = int(ts_score)
        except ValueError:
            return {"Status": "Not yet evaluated"}
        if (hi_score in (0, 1, 2, 3)) or (ts_score in (0, 1, 2, 3)):
            return {"negated": False}
        elif hi_mondo_id == 30 or ts_score == 30:
            return {"negated": False, "mode_of_inheritence": "Autosomal Recessive"}
        elif hi_score == 40 or ts_score == 40:  # another condition for negation
            return {"negated": True}
    else:  # Negate only if no HI disease and TS disease is available for the gene
        return {"negated": True}