import os
import enum
import gzip

from Common.extractor import Extractor
from Common.loader_interface import SourceDataLoader
from Common.node_types import PRIMARY_KNOWLEDGE_SOURCE, NODE_TYPES, SEQUENCE_VARIANT
from Common.prefixes import HGNC  # only an example, use existing curie prefixes or add your own to the prefixes file
from Common.utils import GetData
from datetime import date


# if parsing a tsv or csv type file with columns, use a enum to represent each field
class ClinGenVariantPathogenicityCOLS(enum.IntEnum):
    VARIATION = 0
    CLINVAR_VARIATION_ID = 1
    ALLELE_REGISTRY_ID = 2
    HGVS_EXPRESSIONS = 3
    HGNC_GENE_SYMBOL = 4
    DISEASE = 5
    MONDO_ID = 6
    MODE_OF_INHERITANCE = 7
    ASSERTION = 8
    APPLIED_EVIDENCE_CODES_MET = 9
    APPLIED_EVIDENCE_CODES_NOT_MET = 10
    SUMMARY_OF_INTERPRETATION = 11
    PUBMED_ARTICLES = 12
    EXPERT_PANEL = 13
    GUIDELINE = 14
    APPROVAL_DATE = 15
    PUBLISHED_DATE = 16
    RETRACTED = 17
    EVIDENCE_REPO_LINK = 18
    UUID = 19


##############
# Class: ClinGenVariantPathogenicity  source loader
#
# Desc: Class that loads/parses the ClinGenVariantPathogenicity data.
##############
class ClinGenVariantPathogenicityLoader(SourceDataLoader):
    source_id: str = 'ClinGenVariantPathogenicity'
    # this should be a valid infores curie from the biolink infores catalog
    provenance_id: str = 'infores:clingen'
    # increment parsing_version whenever changes are made to the parser that would result in changes to parsing output
    parsing_version: str = '1.0'
    has_sequence_variants = True  # Flag to use robokop_genetics server to tackle sequence vaiant dataΩ

    def __init__(self, test_mode: bool = False, source_data_dir: str = None):
        """
        :param test_mode - sets the run into test mode
        :param source_data_dir - the specific storage directory to save files in
        """
        super().__init__(test_mode=test_mode, source_data_dir=source_data_dir)

        self.example_url = 'http://erepo.clinicalgenome.org/evrepo/api/classifications/'
        self.example_data_file = 'all?format=tabbed'
        self.data_files = [self.example_data_file]

    def get_latest_source_version(self) -> str:
        # if possible go to the source and retrieve a string that is the latest version of the source data
        latest_version = date.today().strftime("%Y%m%d")
        return latest_version

    def get_data(self) -> bool:
        # get_data is responsible for fetching the files in self.data_files and saving them to self.data_path
        source_data_url = f'{self.example_url}{self.example_data_file}'
        data_puller = GetData()
        print(source_data_url)
        print(self.data_path)
        print(self.example_data_file)
        data_puller.pull_via_http(source_data_url, self.data_path)
        return True

    def parse_data(self) -> dict:
        """
        Parses the data file for graph nodes/edges

        :return: ret_val: load_metadata
        """
        # This is a made up example of how one might extract nodes and edges from a tsv file
        # In this case it's taking the subject ID from column 1 and the object ID from column 3,
        # prepending them with a curie prefix. The predicate comes from column 3. The value in column 4
        # is set as a property on the edge.
        extractor = Extractor(file_writer=self.output_file_writer)
        example_file: str = os.path.join(self.data_path, self.example_data_file)
        print(os.path.join(self.data_path, self.example_data_file))

        with open(example_file, 'rt') as fp:
            extractor.csv_extract(fp,
                                  lambda line: f'CLINVARVARIANT:{line[ClinGenVariantPathogenicityCOLS.CLINVAR_VARIATION_ID.value]}',  # subject id
                                  lambda line: f'{line[ClinGenVariantPathogenicityCOLS.MONDO_ID.value]}',  # object id
                                  lambda line: 'is pathogenic for',  # predicate extractor
                                  lambda line: {NODE_TYPES: SEQUENCE_VARIANT,
                                                'Variation':line[ClinGenVariantPathogenicityCOLS.VARIATION.value],
                                                'HGVS_Gene_Symbol':line[ClinGenVariantPathogenicityCOLS.HGNC_GENE_SYMBOL.value]},  # subject properties
                                  lambda line: {},  # object properties
                                  lambda line: {PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                                                'Assertion' : line[ClinGenVariantPathogenicityCOLS.ASSERTION.value],
                                                'Mode_Of_Inheritance': line[ClinGenVariantPathogenicityCOLS.MODE_OF_INHERITANCE.value],
                                                'Applied_Evidence_Codes_Met':line[ClinGenVariantPathogenicityCOLS.APPLIED_EVIDENCE_CODES_MET.value],
                                                'Applied_Evidence_Codes_Not_Met':line[ClinGenVariantPathogenicityCOLS.APPLIED_EVIDENCE_CODES_NOT_MET.value],
                                                "Summary":line[ClinGenVariantPathogenicityCOLS.SUMMARY_OF_INTERPRETATION.value],
                                                "Pubmed_Articles":line[ClinGenVariantPathogenicityCOLS.PUBMED_ARTICLES.value],
                                                "Expert_Panel":line[ClinGenVariantPathogenicityCOLS.EXPERT_PANEL.value],
                                                "Evidence_Repo_Link":line[ClinGenVariantPathogenicityCOLS.EVIDENCE_REPO_LINK.value],
                                                'Guideline':line[ClinGenVariantPathogenicityCOLS.GUIDELINE.value],
                                                'Approval_Date':line[ClinGenVariantPathogenicityCOLS.APPROVAL_DATE.value],
                                                "Published_Date":line[ClinGenVariantPathogenicityCOLS.APPROVAL_DATE.value]},  # edge properties
                                  comment_character='#',
                                  delim='\t',
                                  has_header_row=True)
        return extractor.load_metadata
