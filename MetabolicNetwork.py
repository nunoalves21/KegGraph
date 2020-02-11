# -*- coding: utf-8 -*-

from MyGraph import MyGraph

class MetabolicNetwork (MyGraph):
    
    def __init__(self, network_type = "metabolite-reaction", split_rev = False):
        MyGraph.__init__(self, {})
        self.net_type = network_type
        self.node_types = {}
        if network_type == "metabolite-reaction":
            self.node_types["metabolite"] = []
            self.node_types["reaction"] = []
        self.split_rev =  split_rev
    
    def add_vertex_type(self, v, nodetype):
        self.add_vertex(v)
        self.node_types[nodetype].append(v)
    
    def get_nodes_type(self, node_type):
        if node_type in self.node_types:
            return self.node_types[node_type]
        else: return None
    
    def load_from_file(self, filename):
        rf = open(filename)
        gmr = MetabolicNetwork("metabolite-reaction")
        for line in rf:
            if ":" in line:
                tokens = line.split(": ")     #alteração de ":" para ": " para não dar conflito com os ids do KEGG
                reac_id = tokens[0].strip()
                gmr.add_vertex_type(reac_id, "reaction")
                rline = tokens[1]
            else: raise Exception("Invalid line:")                
            if "<=>" in rline:
                left, right = rline.split("<=>")
                mets_left = left.split("+")
                for met in mets_left:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    if self.split_rev:
                        gmr.add_vertex_type(reac_id+"_b", "reaction")
                        gmr.add_edge(met_id, reac_id)
                        gmr.add_edge(reac_id+"_b", met_id)
                    else:
                        gmr.add_edge(met_id, reac_id)
                        gmr.add_edge(reac_id, met_id)
                mets_right = right.split("+")
                for met in mets_right:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    if self.split_rev:
                        gmr.add_edge(met_id, reac_id+"_b")
                        gmr.add_edge(reac_id, met_id)
                    else:
                        gmr.add_edge(met_id, reac_id)
                        gmr.add_edge(reac_id, met_id)
            elif "=>" in line:
                left, right = rline.split("=>")
                mets_left = left.split("+")
                for met in mets_left:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    gmr.add_edge(met_id, reac_id)
                mets_right = right.split("+")
                for met in mets_right:
                    met_id = met.strip()
                    if met_id not in gmr.graph:
                        gmr.add_vertex_type(met_id, "metabolite")
                    gmr.add_edge(reac_id, met_id)
            else: raise Exception("Invalid line:")    

        
        if self.net_type == "metabolite-reaction": 
            self.graph = gmr.graph
            self.node_types = gmr.node_types
        elif self.net_type == "metabolite-metabolite":
            self.convert_metabolite_net(gmr)
        elif self.net_type == "reaction-reaction": 
            self.convert_reaction_graph(gmr)
        else: self.graph = {}
        
        
    def convert_metabolite_net(self, gmr):
        for m in gmr.node_types["metabolite"]:
            self.add_vertex(m)
            sucs = gmr.get_successors(m)
            for s in sucs:
                sucs_r = gmr.get_successors(s)
                for s2 in sucs_r:
                    if m != s2:
                        self.add_edge(m, s2)

        
    def convert_reaction_graph(self, gmr): 
        for r in gmr.node_types["reaction"]:
            self.add_vertex(r)
            sucs = gmr.get_successors(r)
            for s in sucs:
                sucs_r = gmr.get_successors(s)
                for s2 in sucs_r:
                    if r != s2: self.add_edge(r, s2)

    ## assessing metabolic potential

    def active_reactions(self, active_metabolites):
        if self.net_type != "metabolite-reaction" or not self.split_rev:
            return None
        res = []
        for v in self.node_types['reaction']:
            preds = set(self.get_predecessors(v))
            if len(preds) > 0 and preds.issubset(set(active_metabolites)):
                res.append(v)
        return res

    def produced_metabolites(self, active_reactions):
        res = []
        for r in active_reactions:
            sucs = self.get_successors(r)
            for s in sucs:
                if s not in res: res.append(s)
        return res

    def all_produced_metabolites(self, initial_metabolites):
        mets = initial_metabolites
        cont = True
        while cont:
            cont = False
            reacs = self.active_reactions(mets)
            new_mets = self.produced_metabolites(reacs)
            for nm in new_mets:
                if nm not in mets:
                    mets.append(nm)
                    cont = True
        return mets

    ## exercise 1 of chapter 14
    def final_metabolites(self):
        res = []
        for v in self.graph.keys():
            if v.startswith('cpd:'):   #alteração de M para cpd devido ao prefixo dos nossos metabolitos
                if len(self.get_predecessors(v)) > 0:
                    if self.get_successors(v) == []:
                        res.append(v)
        return res

    ## exercise 2 of chapter 14
    def shortest_path_product(self, initial_metabolites, target_product):
        if target_product in initial_metabolites:
            return []
        metabs = {}
        for m in initial_metabolites:
            metabs[m] = []
        reacs = self.active_reactions(initial_metabolites)
        cont = True
        while cont:
            cont = False
            for r in reacs:
                sucs = self.get_successors(r)
                preds = self.get_predecessors(r)
                for s in sucs:
                    if s not in metabs:
                        previous = []
                        for p in preds:
                            for rr in metabs[p]:
                                if rr not in previous: previous.append(rr)
                        metabs[s] = previous + [r]
                        if s == target_product:
                            return metabs[s]
                        cont = True
            if cont:
                reacs = self.active_reactions(metabs.keys())
        return None

    def initial_metabolites(self):
        """
        função que vai buscar os metabolitos que não apresenta predecessores
        :return: lista de metabolitos sem predecessores
        """
        res = []
        for v in self.graph.keys():
            if v.startswith('cpd:'):    #cpd devido ao prefixo dos nossos metabolitos
                if len(self.get_successors(v)) > 0:
                    if self.get_predecessors(v) == []:
                        res.append(v)
        return res