import os
import enum
from Common.extractor import Extractor
from Common.loader_interface import SourceDataLoader
from biolink_constants import PRIMARY_KNOWLEDGE_SOURCE, NODE_TYPES, SEQUENCE_VARIANT
from Common.utils import GetData
from Common.utils import LoggingUtil
from datetime import date


# Parsing the columns in the csv file downloaded from clingen source ( https://search.clinicalgenome.org/kb/gene-validity/download )
class ClinGenGeneDiseaseValidityCOLS(enum.IntEnum):
    GENE_SYMBOL = 0
    GENE_ID = 1
    DISEASE_LABEL = 2
    DISEASE_ID = 3
    MOI = 4
    SOP = 5
    CLASSIFICATION = 6
    ONLINE_REPORT = 7
    CLASSIFICATION_DATE = 8
    GCEP = 9


##############
# Class: ClinGenGeneDiseaseValidity  source loader
#
# Desc: Class that loads/parses the ClinGenGeneDiseaseValidity data.
##############
class ClinGenGeneDiseaseValidityLoader(SourceDataLoader):
    source_id: str = "ClinGenGeneDiseaseValidity"
    # this should be a valid infores curie from the biolink infores catalog
    provenance_id: str = "infores:clingen"
    # increment parsing_version whenever changes are made to the parser that would result in changes to parsing output
    parsing_version: str = "1.0"
    # source data
    source_data_url: str = "https://search.clinicalgenome.org/kb/gene-validity/download"

    def __init__(self, test_mode: bool = False, source_data_dir: str = None):
        """
        :param test_mode - sets the run into test mode
        :param source_data_dir - the specific storage directory to save files in
        """
        super().__init__(test_mode=test_mode, source_data_dir=source_data_dir)
        self.data_url = "https://search.clinicalgenome.org/kb/gene-validity/"
        self.gene_disease_data_file = "download"
        self.data_files = [self.gene_disease_data_file]

    def get_latest_source_version(self) -> str:
        # No version is available at the source, using the year_month as versioning proxy
        latest_version = date.today().strftime("%Y%m")
        return latest_version

    def get_data(self) -> bool:
        # get_data is responsible for fetching the files in self.data_files and saving them to self.data_path
        source_data_url = f"{self.data_url}{self.gene_disease_data_file}"
        data_puller = GetData()
        data_puller.pull_via_http(source_data_url, self.data_path)
        return True

    # Created function for normalizing MODE_OF_INHERITANCE of the disease using the moi_lookup dictionary

    moi_lookup = {
        "AD": {
            "ClinGen_label": "AD",
            "normalized_label": "Autosomal Dominant",
            "HPO": "0000006",
        },
        "AR": {
            "ClinGen_label": "AR",
            "normalized_label": "Autosomal Recessive",
            "HPO": "0000007",
        },
        "MT": {
            "ClinGen_label": "MT",
            "normalized_label": "Mitochondrial",
            "HPO": "0001427",
        },
        "SD": {
            "ClinGen_label": "SD",
            "normalized_label": "Semidominant",
            "HPO": "0032113",
        },
        "XL": {"ClinGen_label": "XL", "normalized_label": "X-linked", "HPO": "0001417"},
        "UD": {
            "ClinGen_label": "UD",
            "normalized_label": "Undetermined Mode of Inheritance",
            "HPO": None,
        },
    }

    def moi_normalizer(self, moi, gene, disease):
        try:
            normalized_moi = self.moi_lookup[moi]
        except KeyError:
            normalized_moi = {"label": None, "HPO": None}
            self.logger.info(
                f"No mapping available for {moi} in the moi lookup dictionary for the gene {gene} - disease {disease} pair"
            )
        return normalized_moi

    def parse_data(self) -> dict:
        """
        Parses the data file for graph nodes/edges

        :return: ret_val: load_metadata
        """
        extractor = Extractor(file_writer=self.output_file_writer)
        gene_disease_data_file: str = os.path.join(
            self.data_path, self.gene_disease_data_file
        )

        # Need to incorporate the logic to skip the initial metadata rows
        # either skipped record counter or using loop, rn due the ORION normalization it's working fine
        with open(gene_disease_data_file, "rt") as fp:
            extractor.csv_extract(
                fp,
                lambda line: f"{line[ClinGenGeneDiseaseValidityCOLS.GENE_ID.value]}",  # subject id
                lambda line: f"{line[ClinGenGeneDiseaseValidityCOLS.DISEASE_ID.value]}",  # object id
                lambda line: "gene_associated_with_condition",  # predicate extractor
                lambda line: {},  # subject properties
                lambda line: {},  # object properties
                lambda line: {
                    "MODE_OF_INHERITANCE": self.moi_normalizer(
                        line[ClinGenGeneDiseaseValidityCOLS.MOI.value],
                        line[ClinGenGeneDiseaseValidityCOLS.GENE_ID.value],
                        line[ClinGenGeneDiseaseValidityCOLS.DISEASE_ID.value],
                    ),
                    "CLINGEN_VALIDITY_CLASSIFICATION": line[
                        ClinGenGeneDiseaseValidityCOLS.CLASSIFICATION.value
                    ],
                    "CLINGEN_CLASSIFICATION_DATE": line[
                        ClinGenGeneDiseaseValidityCOLS.CLASSIFICATION_DATE.value
                    ],
                    "CLINGEN_CLASSIFICATION_REPORT": line[
                        ClinGenGeneDiseaseValidityCOLS.ONLINE_REPORT.value
                    ],
                },  # edge properties
                comment_character="#",
                delim=",",
                has_header_row=True,
            )
        return extractor.load_metadata
