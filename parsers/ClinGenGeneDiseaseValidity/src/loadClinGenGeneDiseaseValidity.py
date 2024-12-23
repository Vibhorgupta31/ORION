import os
import enum
import gzip
import re
from Common.extractor import Extractor
from Common.loader_interface import SourceDataLoader
from biolink_constants import PRIMARY_KNOWLEDGE_SOURCE, NODE_TYPES, SEQUENCE_VARIANT
from Common.utils import GetData
from Common.utils import LoggingUtil
import logging
from datetime import date


# Parsing the columns in the csv file downloaded from clingen source
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

    def __init__(self, test_mode: bool = False, source_data_dir: str = None):
        """
        :param test_mode - sets the run into test mode
        :param source_data_dir - the specific storage directory to save files in
        """
        super().__init__(test_mode=test_mode, source_data_dir=source_data_dir)
        # Issue with the file path, need to check this out
        self.data_url = "https://search.clinicalgenome.org/kb/gene-validity/"
        self.gene_disease_data_file = "download"
        self.data_files = [self.gene_disease_data_file]

    def get_latest_source_version(self) -> str:
        # No version is availble at the source, using the date as versioning proxy
        latest_version = date.today().strftime("%Y%m%d")
        return latest_version

    def get_data(self) -> bool:
        # get_data is responsible for fetching the files in self.data_files and saving them to self.data_path
        source_data_url = f"{self.data_url}{self.gene_disease_data_file}"
        data_puller = GetData()
        data_puller.pull_via_http(source_data_url, self.data_path)
        return True

    def parse_data(self) -> dict:
        """
        Parses the data file for graph nodes/edges

        :return: ret_val: load_metadata
        """
        # raw comments **
        # This is a made up example of how one might extract nodes and edges from a tsv file
        # In this case it's taking the subject ID from column 1 and the object ID from column 3,
        # prepending them with a curie prefix. The predicate comes from column 3. The value in column 4
        # is set as a property on the edge.
        extractor = Extractor(file_writer=self.output_file_writer)
        gene_disease_data_file: str = os.path.join(
            self.data_path, self.gene_disease_data_file
        )
        # print(os.path.join(self.data_path, self.gene_disease_data_file)) for_debugging

        # Need to incorporate the logic to skip the intial metadata rows
        # either skipped record counter or using loop, rn the normaliztion version is working fine
        with open(gene_disease_data_file, "rt") as fp:
            extractor.csv_extract(
                fp,
                lambda line: f"{line[ClinGenGeneDiseaseValidityCOLS.GENE_ID.value]}",  # subject id
                lambda line: f"{line[ClinGenGeneDiseaseValidityCOLS.DISEASE_ID.value]}",  # object id
                lambda line: "gene_associated_with_condition",  # predicate extractor
                lambda line: {},  # subject properties
                lambda line: {},  # object properties
                lambda line: {
                    "Mode_of_Inheritance": line[
                        ClinGenGeneDiseaseValidityCOLS.MOI.value
                    ],
                    "Classification": line[
                        ClinGenGeneDiseaseValidityCOLS.CLASSIFICATION.value
                    ],
                    "Classification_Date": line[
                        ClinGenGeneDiseaseValidityCOLS.CLASSIFICATION_DATE.value
                    ],
                    "Classification_Report": line[
                        ClinGenGeneDiseaseValidityCOLS.ONLINE_REPORT.value
                    ],
                },  # edge properties
                comment_character="#",
                delim=",",
                has_header_row=True,
            )
        return extractor.load_metadata


# We can't normalize the Mode of Inheritance or MOI , as the data in the file is not suffice enough to get the MOI
