
__author__ = 'andreap'



class EFOManager(object):

    def __init__(self):
        self.efos = []

    def __contains__(self, item):
        return item.lower() in (n.lower() for n in self.efos)




class Node():

    def __init__(self, id, name, definition = '',synonyms=None):
        self.id = id
        self.name = name
        self.definition = definition
        self.synonyms = synonyms

    def __str__(self):
        return "id:%s | name:%s | definition:%s | synonyms:%s"%(self.id,
                                                   self.name,
                                                   self.definition,
                                                   self.synonyms)




class OBOParser():

    def __init__(self, filename):
        self.filename = filename
        self.efos = []


    def parse(self):
        single_node = []
        return
        '''
        store = False
        for line in file(self.filename):

            if not line.strip():
                store = False
                if single_node:
                    self.efos.append(self._parse_single_node(single_node))
                single_node = []
            if line.startswith('[Term]'):
                store = True
            if store:
                single_node.append(line)
        '''
        
    def _parse_single_node(self, single_node):
        data = dict()
        current_field = ''
        for line in single_node:
            if line.startswith('id'):

                line = line.split(': ')[1]
                synonyms = []
                xrefs = []
                icd9s = []
                data['id'] = line.strip(' \n')
                name = ''
                definition = ''
            elif line.startswith('name:'):

                line = line.split(': ')[1]
                name = line.strip(' \n')
                synonyms.append(name.lower())
            # elif line.startswith('def'):
            #
            #     line = line.split(': ')[1]
            #     definition = line

            elif line.startswith('synonym'):

                line = line.split(': ')[1]
                line = line.split("\"")[1]
                synonyms.append(line.lower())
            elif line.startswith('xref:'):
                line = line.split(': ')[1]
                xrefs.append(line.strip(' \n'))


            elif line.startswith('property_value:'):
                line = line.split(': ')[1]
                line = line.split(' ')
                if line[0] == 'http://www.ebi.ac.uk/efo/ICD9CM_definition_citation' or line[
                    0] == 'http://www.ebi.ac.uk/efo/ICD9_definition_citation':
                    line = line[1].split(':')[1].split('-')

                    for icd9 in line:
                        try:
                            icd9s.append(float(icd9))
                        except ValueError:
                            icd9s.append(icd9)
                elif line[0].endswith('definition_citation'):
                    xrefs.append(line[1])



        data['synonyms'] = synonyms
        data['xref'] = xrefs
        data['icd9'] = icd9s
        return data

    def get_obsolete_efos(self):
        single_node = []
        obsolete_efos = dict()
        store = False
        for line in file(self.filename):

            if not line.strip():
                store = False
                if single_node:
                    new_efo = None
                    obsolete_efo = None
                    for line_node in single_node:
                        if line_node.startswith('id'):
                            line_node = line_node.split(': ')[1]
                            id = line_node.strip(' \n')
                        if line_node.startswith('is_obsolete:'):
                            obsolete_efo = id.replace(':','_')
                        if line_node.startswith('replaced_by:'):
                            try:

                                line_node = line_node.split(': ')[1]
                                line_node = line_node.split('/')
                                new_efo = line_node[4].rstrip()
                            except IndexError as e:
                                new_efo = None

                    if obsolete_efo:
                        obsolete_efos[obsolete_efo] = new_efo


                single_node = []
            if line.startswith('[Term]'):
                store = True
            if store:
                single_node.append(line)
        obsolete_efos['other'] = None
        print(len(obsolete_efos))
        return obsolete_efos