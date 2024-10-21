
"""node, edge, and graph classes and methods for dB graphs"""

from collections import defaultdict
import random


class Node:
    def __init__(self, row, col):
        self.ID = (row, col)
        self.neighbors = []
        self.occupied = False
        self.n_visits = 0

    def add_neighbor(self, neighbor_ID):
        if neighbor_ID not in self.neighbors:
            self.neighbors.append(neighbor_ID)
        self.n_visits += 1

    def remove_neighbor(self, neighbor_ID):
        if neighbor_ID in self.neighbors:
            self.neighbors.remove(neighbor_ID)


class Edge:
    def __init__(self, start, end, weight=1):
        self.start = start
        self.end = end
        self.weight = weight


class DeBruijnGraph():
    """Main class for De Bruijn graphs
    
    Private Attributes:
        graph (defaultdict of lists): Edges for De Bruijn graph
        first_node (str): starting position for traversing the graph
    """

    def __init__(self, input_string, k):
        self.graph = defaultdict(list)  # key corresponds to a node (k-1-mer) 
                                        # value is a list of neighbors (other nodes to which it points)
        self.first_node = ''
        self.nodes = {}     # store nodes by their ID (prefix/suffix)
        self.edges = []     # store edges (kmers)
        self.build_debruijn_graph(input_string, k)
    
    def get_kmers(self):
        """Extract all k-mers from the sequence."""
        for i in range(len(self.sequence) - self.k + 1):
            yield self.sequence[i:i + self.k]

    def add_edge(self, left, right):
        ''' This function adds a new edge to the graph
        
        Args:
            left (str): The k-1 mer for the left edge
            right (str): The k-1 mer for the right edge

        Updates graph attribute to add right to the list named left in defaultdict   
        '''
        edge = Edge(left, right)
        self.edges.append(edge)
        self.nodes[left].add_neighbor(right)
        self.graph[left].append(right)

        return self.graph

    def remove_edge(self, left, right):
        ''' This function removes an edge from the graph
        
        Args:
            left (str): The k-1 mer for the left edge
            right (str): The k-1 mer for the right edge

        Updates graph attribute to remove right from the list named left in defaultdict
        '''
        edge = Edge(left, right)
        if edge in self.edges:
            self.edges.remove(edge)

        if right in self.graph[left]:
            self.graph.remove[right]

        self.nodes[left].remove_neighbor(right)

        return self.graph

        
    def build_debruijn_graph(self, input_string, k):
        """ This function builds a De Buijn graph from a string
        
        Args:
            input_string (str): string to use for building the graph
            k (int): k-mer length for graph construction

        Updates graph attribute to add all valid edges from the string:
            define substring length k and input string
            For each k-length substring of input:
                split k mer into left and right k-1 mer
                add k-1 mers as nodes with a directed edge from left k-1 mer to right k-1 mer"""
        
        """ 1. Extract k-mers from the input sequence.
            2. Use the (k-1)-mer prefixes and suffixes to define nodes.
            3. Create edges between nodes based on overlapping (k-1)-mers.
            4. Track the connections between nodes by using their row/col IDs.
            
            Node class uses (row, col) as ID; we can map each k-1-mer to a unique row/col pair 
            (using indexes or hashing). Each edge will connect nodes whose k-1-mers overlap."""
 
        # Validate input: Ensure k is appropriate for the input string
        if len(input_string) < k:
            raise ValueError("k must be smaller than or equal to the length of the input string.")
        
        # Initialize first node to track the start of the traversal
        first_kmer = input_string[:k]
        self.first_node = first_kmer[:-1]  # Use prefix of first k-mer as the starting point
        
        # extract nodes from k-mers
        for kmer in self.get_kmers(input_string, k):
            prefix = kmer[:-1]
            suffix = kmer[1:]

            # Create nodes if they don't exist
            if prefix not in self.nodes:
                self.nodes[prefix] = Node(prefix)
            if suffix not in self.nodes:
                self.nodes[suffix] = Node(suffix)

            # Create an edge from prefix to suffix
            self.add_edge(prefix, suffix)

        """       
        Example:
        >>> dbg = DeBruijnGraph("this this this is a test", 4)
        >>> print(dbg.graph) #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
        defaultdict(<class 'list'>, {'thi': ['his', 'his', 'his'], 'his': ['is ', 'is ', 'is '], ...)

        """
        
        return self.graph




