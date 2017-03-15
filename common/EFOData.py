
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

    def _parse_single_node(self, single_node):
        data = dict()
        current_field = ''
        for line in single_node:
            if line.startswith('id'):

                line = line.split(': ')[1]
                synonyms = []
                data['id'] = line.strip(' \n')
                name = ''
                definition = ''
            elif line.startswith('name:'):

                line = line.split(': ')[1]
                name = line.strip(' \n')
                synonyms.append(name.lower())
            elif line.startswith('def'):

                line = line.split(': ')[1]
                definition = line

            elif line.startswith('synonym'):

                line = line.split(': ')[1]
                line = line.split("\"")[1]
                synonyms.append(line.lower())

        data['synonyms'] = synonyms
        return data