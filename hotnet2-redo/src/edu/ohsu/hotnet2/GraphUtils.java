package edu.ohsu.hotnet2;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.math3.linear.RealMatrix;

import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.algorithms.filters.FilterUtils;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;

public class GraphUtils{
	private Path currentPath; 
	private String directory; 
	private String fiFile;
	private Path fiFilePath; 
	private String allGenesFile;
	private String allPairsFile;
	private String largestCompFile;
	private String refLargestCompFile;
	private String geneIndexFile;
	private String edgeListFile;
	private FileUtils fu;
	
	public GraphUtils(){
		this.currentPath = Paths.get("");
		this.directory = currentPath.toAbsolutePath().toString();
		this.allGenesFile = "allGenes.txt";
		this.allPairsFile = "allPairs.txt";
		this.largestCompFile = "largestComponentPairs.txt";
		this.refLargestCompFile = "FIsInGene_031516_BigComp.txt";
		this.geneIndexFile = "geneIndexReactome.txt";
		this.edgeListFile = "edgeListReactome.txt";
	}

	public String getFifile(){
		return this.fiFile;
	}
	public void setFifile(String reactomeFileName){
		this.fiFile = reactomeFileName;
	}
		
	/**
	 * Saves 3 files related to the Reactome FI network:
	 * <p>
	 * <ul>
	 * <li>All genes from the Reactome FI network
	 * <li>All interaction pairs from the Reactome FI network
	 * <li>All interaction pairs from the largest component in the Reactome FI network
	 * </ul>
	 * <b>Note:</b> this can also check that the largest component found matches the largest component in a reference file. 
	 * <p>
	 * @param directory	- Directory of Reactome FI network file and place to save files. 	
	 * @param fiFile - Reactome FI network file name.
	 * @param allGenesFile - File name for genes in the Reactome FI network. 
	 * @param allPairsFile - File name for all interaction pairs in the Reactome FI network.
	 * @param largestCompFile - File name for interaction pairs from the largest component.
	 * @throws IOException
	 */
	public void saveReactomeFIGraphFiles(String directory, String fiFile, String allGenesFile, String allPairsFile, String largestCompFile) throws IOException{
		Path fiFilePath = Paths.get(directory, fiFile);
		
		//Get all genes from interaction file and save
		Set<String> allGenes = fu.getAllGenesReactome(fiFilePath);
		fu.saveSetToFile(directory+"/output/", allGenesFile, allGenes);
		
		//Get all interaction pairs from interaction file and save
		Set<String> allPairs = fu.getAllInteractionPairsReactome(fiFilePath);
		fu.saveSetToFile(directory+"/output/", allPairsFile, allPairs);
		
		//Create graph, get largest component, and save the interaction pairs in the largest component
		Graph<String, String> allGenesGraph = createGraph(allGenes, allPairs);
		Graph<String, String> largestComponent = createLargestComponentGraph(allGenesGraph);
		Set<String> genePairs = getPairs(largestComponent);
		fu.saveSetToFile(directory+"/output/", largestCompFile, genePairs);
		
		//Check that interaction pairs in largest component matches reference file
		fu.compareFiles(directory, refLargestCompFile, "/output/"+largestCompFile);
	}
		
	/**
	 * Creates a graph of the whole Reactome FI network.
	 * @param directory - Directory of Reactome FI network file and place to save files.
	 * @param fiFile - Reactome FI network file name.
	 * @return a graph of the whole Reactome FI network.
	 * @throws IOException
	 */
	public Graph<String, String> createReactomeFIGraph(String directory, String fiFile) throws IOException{
//		FileUtils fu = new FileUtils();
		Path fiFilePath = Paths.get(directory, fiFile);
		Set<String> allGenes = fu.getAllGenesReactome(fiFilePath);
		Set<String> allPairs = fu.getAllInteractionPairsReactome(fiFilePath);
		Graph<String, String> allGenesGraph = createGraph(allGenes, allPairs);
		return allGenesGraph; 
	}
	
	/**
	 * Creates a graph.
	 * @param genes - Set of genes in the graph.
	 * @param pairs - Set of interaction pairs in the graph.
	 * @return a undirected graph. 
	 */
	public Graph<String,String> createGraph (Set<String> genes, Set<String> pairs){
		Graph<String, String> graph = new UndirectedSparseGraph<String, String>();	
		for (String g: genes)
			graph.addVertex(g);
		if(!pairs.isEmpty()){
			for (String p: pairs){ 
				String[] row = p.split("\t");
				String gene1 = row[0];
				String gene2 = row[1];
				graph.addEdge(p ,gene1, gene2); 
			}			
		}
		return graph;
	}

	/**
	 * Creates a directed graph.
	 * @param genes - Set of genes in the graph.
	 * @param pairs - Set of interaction pairs in the graph.
	 * @return a directed graph. 
	 */
	private Graph<String,String> createGraphDirected (Set<String> genes, Set<String> pairs){
		Graph<String, String> graph = new DirectedSparseGraph<String, String>();	
		for (String g: genes)
			graph.addVertex(g);
		if(!pairs.isEmpty()){
			for (String p: pairs){ 
				String[] row = p.split("\t");
				String gene1 = row[0];
				String gene2 = row[1];
				graph.addEdge(p ,gene1, gene2); 
			}			
		}
		return graph;
	}
	
	/**
	 * Gets the largest component in a graph.
	 * @param graph - Graph used to find the largest component.
	 * @return a ordered set of genes in the largest component.
	 */
	private Set<String> getLargestComponent(Graph<String, String> graph){
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		Set<String> largestComponent = new HashSet<String>();
		for (Set<String> c: components){	
			if(c.size() > largestComponent.size())
				largestComponent = c; 
		}
		return largestComponent; 
	}
		
	/**
	 * Creates a graph of only the largest component.
	 * @param graph - Graph of network.
	 * @return	a graph of only the largest component within provided network.
	 */
	public Graph<String, String> createLargestComponentGraph (Graph<String, String> graph){
		Set<String> genes = getLargestComponent(graph);
		Graph<String, String> largestComponentGraph = FilterUtils.createInducedSubgraph(genes,graph);	
		return largestComponentGraph; 
	}
	
	/**
	 * Gets interaction pairs from a graph.
	 * @param graph - graph for extracting edges.
	 * @return a set containing the graph's interaction pairs.  
	 */
	public Set<String> getPairs(Graph<String, String> graph){
		Set<String> pairs = new HashSet<String>();
		for (String v: graph.getVertices()){
			for (String e: graph.getEdges()){
				//note: graph.Endpoints(e) stores edges as Pair<L,R> (ex. <gene1, gene2>), so the first "<", ", ", and the last ">" must be removed
				String tempString = graph.getEndpoints(e).toString().replaceFirst("<", "").replace(", ", "\t").replaceFirst("[>]$","");	
				pairs.add(tempString);
			}
			break;
		}
		return pairs; 
	}
	
	/**
	 * Gets an ordered set of genes in the graph.
	 * @param graph	- Graph used to obtain genes from. 
	 * @return a ordered set of genes from the graph.
	 */
	public SortedSet<String> getGeneGraphSet (Graph<String, String> graph){
		SortedSet<String> geneSet = new TreeSet<String>();
		for (String vertex: graph.getVertices())
			geneSet.add(vertex);
		return geneSet;
	}
	
	/**
	 * Gets a list of degree for each gene in the provided set.  
	 * @param graph	- Graph used to obtain degree from for each gene.
	 * @param geneSet - Set of genes from the graph used to determine the list order of gene degree.  
	 * @return a list containing the degree for each gene in the provided set.
	 */
	public List<Integer> getGeneDegreeList(Graph<String,String> graph, Set<String> geneSet){
		List<Integer> degreeList = new ArrayList<Integer>();
		List<String> geneList = new ArrayList<String>(geneSet);
		for (String g: geneList){
			degreeList.add(graph.getOutEdges(g).size());
		}
		return degreeList;
	}
	
	/**
	 * Converts the provided RealMatrix into a directed graph. 
	 * @param matrix - Matrix to convert into graph.
	 * @param geneSet - Set of genes used to determine matrix order.
	 * @return a directed graph generated from the matrix
	 */
	public Graph<String, String> covertMatrixToDirectedGraph (RealMatrix matrix, SortedSet<String> geneSet){
		List<String> geneList = new ArrayList<String>(geneSet);
		Set<String> genes = new TreeSet<String>();
		Set<String> pairs = new TreeSet<String>();
		int dim = matrix.getRowDimension();
		for (int i=0; i<dim; i++){
			for (int j=0; j<dim; j++){
				double value = matrix.getEntry(i, j);
				if (value > 0){
					String gene1 = Integer.toString(i);
					String gene2 = Integer.toString(j);
					genes.add(gene1);
					genes.add(gene2);
					pairs.add(gene1 + "\t" + gene2);
				}
			}
		}	
		Graph<String, String> graph = createGraphDirected(genes, pairs); 
		return graph; 
	}
	
	/**
	 * Gets a degree hash map of a graph using gene as key and degree as value.
	 * <p>
	 * Based on Python NetworkX's <a href="https://networkx.readthedocs.io/en/stable/_modules/networkx/classes/graph.html#Graph.degree">Graph.degree()</a>.
	 * @param graph	- Graph used to create degree HashMap.
	 * @return a HashMap with gene as key and gene's degree as value.
	 */
	public HashMap<String, Integer> getDegreeMap(Graph<String,String> graph){
		HashMap<String, Integer> degreeMap = new HashMap<String, Integer>();
		for (String v: graph.getVertices())
			degreeMap.put(v, graph.getOutEdges(v).size());
		return degreeMap;
	}
	
}