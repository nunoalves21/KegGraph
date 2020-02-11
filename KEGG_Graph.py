
import urllib.request
from Bio.KEGG.KGML import KGML_parser
from Bio.KEGG import REST, Compound
from MetabolicNetwork import MetabolicNetwork

class KEGG_Graph(MetabolicNetwork):

    def __init__(self, pathway_id, network_type = "metabolite-reaction", split_rev = False):
        MetabolicNetwork.__init__(self, network_type, split_rev)
        self.reactions = {}     #dicionário com os id's do KEGG das reações como chaves e respetivos nomes como valores
        self.metabolites = {}       #igual ao de cima mas para os metabolitos
        self.name = pathway_id
        self.getKGML()
        self.xml_to_txt()
        self.load_from_file(self.name + '.txt')

    def getKGML(self):
        """
        função que obtem ficheiro xml da pathway através do seu id
        """
        url = 'http://rest.kegg.jp/get/' + self.name + '/kgml'      #url da pathway
        file = self.name + '.xml'
        urllib.request.urlretrieve(url, file)

    def xml_to_txt(self):
        """
        função que converte o ficheiro xml gerado com a função self.getKGML() para a conformação das reações dada nas aulas
        """
        pathway = KGML_parser.read(open(self.name + '.xml', 'r'))
        reactions = pathway._reactions
        f = open(self.name + '.txt', 'w+')
        for r_id in reactions.keys():
            reaction_names = reactions[r_id]._names
            substrates_ids_t = pathway._reactions[r_id]._substrates
            products_ids_t = pathway._reactions[r_id]._products
            while reaction_names != []:         #para o mesmo id por vezes existe mais do que um nome de reação, por isso, este while serve para percorrer todos os nomes existentes
                substrates_ids = list(substrates_ids_t)
                products_ids = list(products_ids_t)
                name = reaction_names.pop()
                self.dictionary(name)       #para adicionar o nome da reação ao dicionário criado
                line = name + ': '
                while substrates_ids:       #mesmo fundamento do while anterior mas para os substratos
                    sub = substrates_ids.pop()
                    self.dictionary(pathway.entries[sub]._names[0])
                    if len(substrates_ids)>0:
                        line += str(pathway.entries[sub]._names[0]) + ' + '
                    else:
                        line += str(pathway.entries[sub]._names[0]) + ' '
                if pathway._reactions[r_id].type == "reversible":   #adicionar os simbolos definidos para reversivel / irreversivel
                    line += '<=> '
                else:
                    line += '=> '
                while products_ids:     #mesmo fundamento dos whiles anteriores mas para os produtos
                    prod = products_ids.pop()
                    self.dictionary(pathway.entries[prod]._names[0])
                    if len(products_ids)>0:
                        line += str(pathway.entries[prod]._names[0]) + ' + '
                    else:
                        line += str(pathway.entries[prod]._names[0]) + '\n'
                f.write(line)
        f.close()


    def dictionary(self, name):
        """
        função que cria um dicionário com os id's do KEGG das reações/metabolitos como chaves e respetivos nomes como valores
        :param name: id do KEGG para a reação/substrato
        """
        if name not in self.metabolites.keys() and name not in self.reactions.keys():
            comp = Compound.parse(REST.kegg_get(name))
            for c in comp:
                names = []
                for n in c.name:
                    names.append(n.lower())
                if name.startswith('cpd:'):
                    self.metabolites[name] = names
                else:
                    self.reactions[name] = names

    def find_metabolite(self, name):
        """
        função que retorna o id do KEGG a partir do respetivo nome do metabolito
        :param name: nome do metabolito
        :return k: id KEGG do metabolito
        """
        for k in self.metabolites.keys():
            for v in self.metabolites[k]:
                if v == name.lower():
                    return k
        return "Metabolite not found"

    def find_reaction(self, name):
        """
        similar à find_metabolite mas para reações
        :param name: nome da reação
        :return k: id KEGG da reação
        """
        for k in self.reactions.keys():
            for v in self.reactions[k]:
                if v == name.lower():
                    return k
        return "Reaction not found"

    def find_idname(self, id):
        """
        função que a partir do id do KEGG retorna o nome do respetivo metabolito/ reação
        :param id: id do KEGG
        :return id, self.metabolites[km][0], self.reactions[kr}[0]: id do KEGG, nome por extenso do metabolito,
        nome por extenso da reação
        """
        if id.endswith('_b'):
            id = id[:-len('-b')]
        for km in self.metabolites.keys():
            if km == id:
                if self.metabolites[km] != []:
                    return self.metabolites[km][0]  #retorna o primeiro nome por extenso que encontra
                else: return id            #caso não encontre nome, retorna novamente o id do KEGG
        for kr in self.reactions.keys():
            if kr == id:
                if self.reactions[kr] != []:
                    return self.reactions[kr][0]
                else: return id
        return id

    def rename(self, result):
        """
        função que percorre um resultado(lista, lista de tuplos, valor, lista de listas, etc.) e substitui cada id do KEGG pelo respetivo nome (caso o encontre)
        :param result: resultado obtido em outras funções
        :return: mesmo resultado mas com os nomes substituídos
        """
        if isinstance(result, str):
            return self.find_idname(result)
        else:
            for i in range(len(result)):
                if isinstance(result[i], str):
                    result=list(result)
                    result[i] = self.find_idname(result[i])
                else:
                    for j in range(len(result[i])):
                        if isinstance(result[i][j], str):
                            result[i] = list(result[i])
                            result[i][j] = self.find_idname(result[i][j])
        return result






if __name__ == "__main__":
    kegg=KEGG_Graph('hsa00010', split_rev=True)
    kegg.print_graph()

    e=kegg.reachable_with_dist('cpd:C00036')
    print(kegg.rename(e))
    f=kegg.shortest_path('cpd:C00036','cpd:C00221')
    print(kegg.rename(f))
    print(kegg.find_idname('rn:R03270'))




