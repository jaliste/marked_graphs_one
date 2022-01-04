from collections import defaultdict

def core_of_graph(T):
    if type(T) == WeightedGraph:
        T = T.copy()
    else:
        T = WeightedGraph(T.copy())
    leaves = [v for v in T.vertices(sort=False) if (T.degree(v)==1 and T.vertex_weights[v]==(1,0))]
    # If we have only one edge 1-1, we can only contract one vertex
    if T.num_verts()==2 and len(leaves)==2:
        leaves = leaves[:-1]
    internal = set(T.vertices(sort=False)).difference(set(leaves))
    for v in internal:
        w,d = T.vertex_weights[v]
        T.vertex_weights[v] = (len(set(T.neighbors(v)).intersection(leaves)) + w,d)
    T.delete_vertices(leaves)
    return T

class WeightedGraph (Graph):
    def __init__(self, data=None, pos=None, loops=True, format=None,
                 weighted=None, implementation='c_graph',
                 data_structure="sparse", vertex_weights=None, add_function = None,
                 vertex_labels=True, name=None, multiedges=True,
                 convert_empty_dict_labels_to_None=None, sparse=True,
                 immutable=False):
        Graph.__init__(self, data, pos, loops, format, weighted,
                       data_structure, vertex_labels, name,multiedges,
                       convert_empty_dict_labels_to_None, sparse, immutable)
        if vertex_weights is None:
            vertex_weights = defaultdict(lambda: (1,0))
        if add_function is None:
            add_function = lambda x,y:(x[0]+y[0],x[1]+y[1]+1)
        self.vertex_weights = vertex_weights
        self.add_function = add_function
        self.vertex_map = dict()
        for v in self.vertices(sort=False):
            self.vertex_map[v] = Set({v})
    def is_proper(self):
        for v in self.vertices(sort=False):
            if self.vertex_weights[v]==(1,0):
                return False
        return True
    def weights_multiset(self):
        return [self.vertex_weights[v] for v in self.vertices(sort=False)]
 
    def label_good_leaves(self):
        
        if self.num_edges()==1:
            if self.is_proper():
                self.leaves0 = self.edges(sort=False)
                
            else:
                self.leaves0 = []
            return
                
        self.leaves0 = self.edges_incident([v for v in self.vertices(sort=False) if self.degree(v)==1])
        
        non_leaf_vertices = [v for v in self.vertices(sort=False) if self.degree(v)>1]
        
        if len(non_leaf_vertices)==1 and self.vertex_weights[non_leaf_vertices[0]]==(1,0):
            #It is a star with a (1,0) vertex. Last leaf cannot be removed!
            self.leaves0 = self.leaves0[:-1]
            #for l in leaves:
            #    self.set_edge_label(l[0],l[1],i)
            #    i+=1
            #    work.delete_vertex(l[0])
        
    def label_edges(self):
        work = self.copy()
        i = 1
        self.num_leaves = []
        self.leaves = []
        while work.num_verts()>1:
            if work.num_edges()>1:
                leaves = [(v,work.neighbors(v)[0]) for v in work.vertices(sort=False) if work.degree(v)==1]
            else:
                leaves = work.edges(sort=False)
            
            self.num_leaves.append(len(leaves))
            self.leaves.append(leaves)
            for l in leaves:
                self.set_edge_label(l[0],l[1],i)
                i+=1
                work.delete_vertex(l[0])
    def dissolve_contractible_leaves(self):
        return self.contract_edges(self.contractible_leaves(),add_func=undot_sum,copy=False)

            
    def copy(self,immutable=False):
        T = WeightedGraph(Graph.copy(self),vertex_weights=self.vertex_weights.copy(),immutable=immutable)
        return T
    def plot2(self):
        from sage.misc.element_with_label import ElementWithLabel
        return self.relabel(lambda v: ElementWithLabel(v,self.vertex_weights[v]), inplace = False).plot(edge_labels=True,vertex_size=550)
        
        
    
    def d_part(self):
        weights = [self.vertex_weights[v] for v in self.vertices(sort=False)]
        return tuple(sorted(weights,reverse=True))
    def __str__(self):
        return str(self.d_part())
    def __repr__(self):
        return str(self.d_part())
    def __unicode__(self):
        return str(self.d_part())
    def map_vertex(self,v):
        try:
            return list(v)
        except TypeError:
            return [v]
    def total_dots(self):
        suma = 0
        for v in self.vertices(sort=False):
            suma += self.vertex_weights[v][1]
        return suma + self.num_edges()
    def total_weight(self):
        suma = 0
        for v in self.vertices(sort=False):
            suma += self.vertex_weights[v][0]
        return suma
    def contractible_leaves(self):
        if self.num_edges()==1:
            if self.is_proper():
                return []
            else:
                return self.edges(sort=False)
        return [(v,self.neighbors(v)[0]) for v in self.vertices(sort=False) if self.degree(v)==1 and self.vertex_weights[v]==(1,0)]
    
    def contract_edges2(self, edges, add_func=None, copy=True):
        if add_func is None:
            add_func = self.add_function
        if copy:
            H = self.copy()
        else:
            H = self

        if len(set(len(e) for e in edges)) > 1:
            raise ValueError("edge tuples in input should have the same length")
        edge_list = []
        vertices = set()
        for e in edges:
            # try to get the vertices and label of e as distinct variables
            u = e[0]
            v = e[1]
            label = H.edge_label(u,v)
            
            if H.has_edge(u, v):
                edge_list.append((u, v))
                vertices.add(u)
                vertices.add(v)
        if not edge_list:
            return H

        # implementation of union_find using DisjointSet
        from sage.sets.disjoint_set import DisjointSet
        DS = DisjointSet(H.vertex_iterator())
        
        for u, v in edge_list:
            DS.union(u, v)
        H.delete_edges(edge_list)

        edges_incident = []
        vertices = [v for v in vertices if v != DS.find(v)]
        
        for v in vertices:
            edges_incident.extend(H.edges_incident(v, sort=False))
            r = DS.find(v)
            H.vertex_weights[r] = add_func(H.vertex_weights[r],H.vertex_weights[v])
            H.delete_vertex(v)

        for (u, v,label) in edges_incident:
            root_u = DS.find(u)
            root_v = DS.find(v)
            
            if root_v != root_u or H.allows_loops():
                H.add_edge(root_u, root_v)
        return H
    
    
    def contract_edges (self, list_ee, add_func=None, copy=True):
        if add_func is None:
            add_func = self.add_function
        if copy:
            H = self.copy()
        else:
            H = self
        for e in list_ee:
            u0 = e[0]
            v0 = e[1]

            count = 0
            for v in H.vertices(sort=False):
                vv = H.map_vertex(v)
                if u0 in vv:
                    u0 = v
                    count +=1
                if v0 in vv:
                    v0 = v
                    count +=1
                if count == 2:
                    break
            u = (H.map_vertex(u0))
            v = (H.map_vertex(v0))
            new_v = tuple(u+v)
            tu = tuple(u)
            tv = tuple(v)
            #print "Contracting",tu,tv," into ", new_v
            H.add_vertex(new_v)
            H.delete_edge(u0,v0)
            H.add_edges([(n,new_v,H.edge_label(n,u0)) for n in H.neighbors(u0)])
            H.add_edges([(n,new_v,H.edge_label(n,v0)) for n in H.neighbors(v0)])
            H.delete_vertices([u0,v0])
            H.vertex_weights[new_v] = add_func(H.vertex_weights[u0], H.vertex_weights[v0])

            del H.vertex_weights[u0]

            del H.vertex_weights[v0]
            
           
        return H
    def Dpolynomial(self):
        Edges = self.edges(sort=False)
        d = len(Edges)
        N = self.total_weight()
       
        self.Lattice = posets.BooleanLattice(len(Edges))
        self.LatEdges = { e: getEdgeSetFromLattice(Edges,e) for e in self.Lattice }
        self.P = { e: self.contract_edges2(self.LatEdges[e]).copy(immutable=True) for e in self.Lattice}
        gens = getGenerators(N,d)
        R = PolynomialRing(ZZ, [_tuple2gen(*i) for i in gens])
        X = R(0)
        for x in self.Lattice:
            d_part = self.P[x]
            X +=NStd2(R,*zip(*d_part.d_part()))
            
        return X

def getEdgeSetFromLattice(E,number):
    return [ E[i] for i,v in enumerate(number.digits(2)) if v==1]
def N_pol(F,a,dots):
    return sum([(-1)**i*int(binomial(dots,i))*F('x1_0')**i*F('x'+str(a-i)+'_0') for i in range(dots+1)])
def NStd(F,a):
    res = F(1)
    for i in range(len(a)):
        try:
            L = a[i][0] 
            k = a[i][1]
        except TypeError:
            L = a[i]
            k = 0
        res *= N_pol(F,L,k)
    return res
def NStd2(F, t, p):
    aux = []
    for i in range(len(t)):
        aux.append((t[i],p[i]))
    return NStd(F, aux)



