import os
import enum
import gzip
from Common.extractor import Extractor
from Common.loader_interface import SourceDataLoader
from biolink_constants import PRIMARY_KNOWLEDGE_SOURCE, NODE_TYPES, SEQUENCE_VARIANT
from Common.utils import GetData
from datetime import date


# Parsing the columns in the tsv file downloaded from clingen source
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
# Class: ClinGenDosageSensitivity  data loader
# Desc: Class that loads/parses the ClinGen Dosage Sensitivity  data.
###########
class ClinGenDosageSensitivityLoader(SourceDataLoader):
    source_id: str = "ClinGenDosageSensitivity"
    provenance_id: str = "infores:clingen"  # From BioLink infores catalog
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
        # No version is availble at the source, using the year_month when the code was run as versioning proxy
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
        extractor = Extractor(file_writer=self.output_file_writer)
        dosage_sensitivity_gene_file: str = os.path.join(
            self.data_path, self.cligen_dosage_sensitivity_gene_file
        )
        dosage_sensitivity_region_file: str = os.path.join(
            self.data_path, self.clingen_dosage_sensitivity_region_file
        )
        # Gene file processing
        ## Haploinsufficiency edges parsing
        with open(dosage_sensitivity_gene_file, "rt") as fp:
            extractor.csv_extract(
                fp,
                lambda line: f"NCBIGene:{line[ClinGenDosageSensitivityCOLS.GENE.value]}",
                # subject id
                lambda line: (
                    "MONDO:0700096"
                    if line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value] == ""
                    else line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value]
                ),  # object id
                lambda line: "gene associated with condition",  # predicate extractor
                lambda line: {},  # subject properties
                lambda line: {},  # object properties
                lambda line: dict(
                    {
                        PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                        "Haploinsufficiency_Description": line[
                            ClinGenDosageSensitivityCOLS.HI_DESCRIPTION.value
                        ],
                        "Haploinsufficiency_Score": line[
                            ClinGenDosageSensitivityCOLS.HI_SCORE.value
                        ],
                    },
                    **get_edge_properties(
                        line[ClinGenDosageSensitivityCOLS.HI_SCORE.value],
                        line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value],
                    ),
                ),
                # edge properties
                comment_character="#",
                delim="\t",
                has_header_row=True,
            )
        ## Triplosensitivity edges parsing
        with open(dosage_sensitivity_gene_file, "rt") as fp:
            extractor.csv_extract(
                fp,
                lambda line: f"NCBIGene:{line[ClinGenDosageSensitivityCOLS.GENE.value]}",
                # subject id
                lambda line: (
                    "MONDO:0700096"
                    if line[ClinGenDosageSensitivityCOLS.TS_DISEASE.value] == ""
                    else line[ClinGenDosageSensitivityCOLS.TS_DISEASE.value]
                ),  # object id
                lambda line: "gene associated with condition",  # predicate extractor
                lambda line: {},  # subject properties
                lambda line: {},  # object properties
                lambda line: dict(
                    {
                        PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                        "Triplosensitivity_Description": line[
                            ClinGenDosageSensitivityCOLS.TS_DESCRIPTION.value
                        ],
                        "Triplosensitivity_Score": line[
                            ClinGenDosageSensitivityCOLS.TS_SCORE.value
                        ],
                    },
                    **get_edge_properties(
                        line[ClinGenDosageSensitivityCOLS.TS_SCORE.value],
                        line[ClinGenDosageSensitivityCOLS.TS_DISEASE.value],
                    ),
                ),
                # edge properties
                comment_character="#",
                delim="\t",
                has_header_row=True,
            )
        # Region file processing
        ## Haploinsufficiency edges parsing
        with open(dosage_sensitivity_region_file, "rt") as fp:
            extractor.csv_extract(
                fp,
                lambda line: f"NCBIGene:{line[ClinGenDosageSensitivityCOLS.REGION.value]}",
                # subject id
                lambda line: (
                    "MONDO:0700096"
                    if line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value] == ""
                    else line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value]
                ),  # object id
                lambda line: "region associated with condition",  # predicate extractor
                lambda line: {},  # subject properties
                lambda line: {},  # object properties
                lambda line: dict(
                    {
                        PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                        "Haploinsufficiency_Description": line[
                            ClinGenDosageSensitivityCOLS.HI_DESCRIPTION.value
                        ],
                        "Haploinsufficiency_Score": line[
                            ClinGenDosageSensitivityCOLS.HI_SCORE.value
                        ],
                    },
                    **get_edge_properties(
                        line[ClinGenDosageSensitivityCOLS.HI_SCORE.value],
                        line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value],
                    ),
                ),
                # edge properties
                comment_character="#",
                delim="\t",
                has_header_row=True,
            )
        ## Triplosensitivity edges parsing
        with open(dosage_sensitivity_region_file, "rt") as fp:
            extractor.csv_extract(
                fp,
                lambda line: f"NCBIGene:{line[ClinGenDosageSensitivityCOLS.REGION.value]}",
                # subject id
                lambda line: (
                    "MONDO:0700096"
                    if line[ClinGenDosageSensitivityCOLS.TS_DISEASE.value] == ""
                    else line[ClinGenDosageSensitivityCOLS.TS_DISEASE.value]
                ),  # object id
                lambda line: "region associated with condition",  # predicate extractor
                lambda line: {},  # subject properties
                lambda line: {},  # object properties
                lambda line: dict(
                    {
                        PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                        "Triplosensitivity_Description": line[
                            ClinGenDosageSensitivityCOLS.TS_DESCRIPTION.value
                        ],
                        "Triplosensitivity_Score": line[
                            ClinGenDosageSensitivityCOLS.TS_SCORE.value
                        ],
                    },
                    **get_edge_properties(
                        line[ClinGenDosageSensitivityCOLS.TS_SCORE.value],
                        line[ClinGenDosageSensitivityCOLS.TS_DISEASE.value],
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
# Scoring values are from https://clinicalgenome.org/site/assets/files/6428/dosage_sop-scoring-1.pdf
# On brief score stands for : 0: No evidence ; 1: Litle evidence ; 2: Emerging evidence ; 3: Sufficient evidence ; 30: gene / regiion associated with AR condition ; 40: Evidence against or dosage sensitivity unlikely ; -1 : Not planned to be evaluated


def get_edge_properties(score, disease_id):
    if disease_id != "":
        try:
            score = int(score)
        except ValueError:
            return {"Status": "Not yet evaluated"}
        if score in (0, 1, 2, 3):
            return {"negated": False}
        elif score == 30:
            return {"negated": False}
        elif score == 40:  #  condition for negation
            return {"negated": True}
    else:  # Negate only if no HI disease and TS disease is available for the gene
        return {"negated": True}
