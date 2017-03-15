import csv
from common.EFOData import OBOParser
from common.HGNCParser import GeneParser
import json
from settings import Config

class Phewas(object):

    def __init__(self, phewas_id, snp = None, ensg_id = None, efo_id = None, phenotype = None, gene_name = None):
        self.ensg_id = ensg_id
        self.efo_id = efo_id
        self.snp = snp
        self.phenotype = phenotype
        self.gene_name = gene_name
        self.phewas_id = phewas_id

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__,
                          sort_keys=True, indent=4)


class CSVParser(object):

    def __init__(self, filename):
        self.filename = filename


    def parse(self):
        with open(self.filename, "rb") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                yield row


class PhewasProcessor(object):
    def __init__(self):
        self.genes = dict()
        self.efos = list()





    def find_gene(self,gene_name):

        #matched_gene = filter(lambda gene:gene.approved_symbol == gene_name,self.genes )
        matched_gene = self.genes.get(gene_name)
        # for gene in  matched_gene:
        #     print 'Gene -->' + gene.ensembl_gene_id
        # if matched_gene:
        #     return matched_gene[0].ensembl_gene_id


    def find_disease(self,phewas_phenotype):

        #matched_disease = filter(lambda disease: disease.name.lower() == phewas_phenotype.lower() or phewas_phenotype.lower() in (n.lower() for n in disease.synonyms), self.diseases)
        matched_disease = [efo_dict['id'] for efo_dict in self.efos if phewas_phenotype.lower() in efo_dict['synonyms'] ]

        # for phenotype in matched_disease:
        #     print 'Disease --> {}'.format(phenotype)
        if matched_disease:
            return matched_disease[0]


    def convert_phewas_catalog_evidence_json(self):

        obo_parser = OBOParser('../efo.obo')
        obo_parser.parse()
        self.efos = obo_parser.efos


        gene_parser = GeneParser()
        gene_parser._get_hgnc_data_from_json()
        self.genes = gene_parser.genes

        csv_parser = CSVParser('../phewas-catalog.csv')
        counter = 1
        my_list = []
        missing_efo_fieldnames = ['phenotype', 'similar_efo']
        fieldnames = ['phenotype', 'efo_id', 'gene_name','ensg_id','snp']
        with open('../missing_efo.csv', 'w') as out_missing_csv , open('../phewas_efo_ensg.csv', 'w') as out_csv:
            writer = csv.DictWriter(out_csv, fieldnames)
            writer.writeheader()
            missing_efo_writer = csv.DictWriter(out_missing_csv, missing_efo_fieldnames)
            missing_efo_writer.writeheader()
            for phewas_row in csv_parser.parse():
                phewas_obj = Phewas(counter)
                # print '------------------------------------------------------------------------------'
                # print 'Phewas --> ' + phewas_row['phewas phenotype']


                phewas_obj.ensg_id =  self.genes.get(phewas_row['gene_name'])#self.find_gene(phewas_row['gene_name'])
                phewas_obj.efo_id = self.find_disease(phewas_row['phewas phenotype'])
                phewas_obj.snp = phewas_row['snp']
                phewas_obj.phenotype = phewas_row['phewas phenotype']
                phewas_obj.gene_name = phewas_row['gene_name']

                if phewas_obj.efo_id :
                    inner_dict = dict(zip(fieldnames, [phewas_obj.phenotype, phewas_obj.efo_id,phewas_obj.gene_name,phewas_obj.ensg_id,phewas_obj.snp]))
                    writer.writerow(inner_dict)
                elif phewas_obj.efo_id is None:
                    inner_dict = dict(zip(missing_efo_fieldnames, [phewas_row['phewas phenotype'], '']))
                    missing_efo_writer.writerow(inner_dict)

                # batch the write process
                # profile pycharm

                #json.dump(phewas_obj.toJSON(),out_json)



def remove_dup():


    with open('../missing_efo.csv', 'r') as in_file, open('../missing_efos.csv', 'w') as out_file:
        seen = set()  # set for fast O(1) amortized lookup
        reader = csv.reader(in_file, dialect=csv.excel_tab)
        for row in reader:
            if row[0] in seen: continue  # skip duplicate

            seen.add(row[0])
            out_file.write(row[0])
            out_file.write('\n')

def main():
    phewas_processor = PhewasProcessor()
    phewas_processor.convert_phewas_catalog_evidence_json()
    #remove_dup()
    test = 'Hypoglycemia'


if __name__ == "__main__":
    main()
