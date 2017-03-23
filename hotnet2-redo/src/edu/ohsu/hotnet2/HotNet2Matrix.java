package edu.ohsu.hotnet2;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.DefaultRealMatrixPreservingVisitor;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SparseRealMatrix;
import org.jgrapht.alg.StrongConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraphModified;
import org.jgrapht.graph.DefaultEdge;
import org.junit.Test;
import org.ojalgo.matrix.PrimitiveMatrix;

import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.DirectedGraph;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.Graph;

public class HotNet2Matrix {
	
	public HotNet2Matrix(){
	
	}
	
	@Test
	public void testPrototype2HotNet2Algo() throws IOException, InterruptedException{
		Path currentPath = Paths.get("");
		String directory = currentPath.toAbsolutePath().toString();
		double beta = 0.25;
		double delta = 0.001;
		int numPermutation = 100;
		String dir =  directory+"/prototype2/";
		String indexFile = "prototype2_index_java";
		String edgeListFile = "prototype2_edgelist_java";
		String heatScoreFile = "prototype2_heatScore_java";
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();
		//Create set of genes in network and HashMap of gene index to allow gene names in graph instead of numbers
		SortedSet<String> genes = fu.getAllGenesPY(dir, indexFile, "\t");
		HashMap<String, String> geneIndexMap = fu.getGeneIndexNamePY(dir, indexFile, "\t");
		//Create graph from file and check there is only 1 component
		Set<String> pairs = fu.getAllInteractionPairsPY(dir, edgeListFile, geneIndexMap, false, "\t");
		Graph<String, String> graph = gu.createGraph(genes, pairs);
		Set<Set<String>> components = wcc.transform(graph);
		if (components.size() != 1){
			//throw new IllegalArgumentException("Provided permuted graph " + edgeListFile + " is not connected, it has " + components.size() + " components.");
			System.out.println("Provided permuted graph " + edgeListFile + " is not connected, it has " + components.size() + " components.");
		}
		
		SortedSet<String> geneSet = gu.getGeneGraphSet(graph);
				
		java.util.concurrent.TimeUnit.SECONDS.sleep(10);
		
		//Comment out the functions not being timed
/*		long start = System.currentTimeMillis();
		OOA(graph, geneSet, beta, delta, dir, heatScoreFile);
		long end = System.currentTimeMillis();
		System.out.println("total time\t OOA: " + (end-start));
*/				
/*		long start = System.currentTimeMillis();
		OOO(graph, geneSet, beta, delta, dir, heatScoreFile);
		long end = System.currentTimeMillis();
		System.out.println("total time\t OOO: " + (end-start));
*/		
/* start = System.currentTimeMillis();
		OAA(graph, geneSet, beta, delta, dir, heatScoreFile);
		long end = System.currentTimeMillis();
		System.out.println("total time\t OAA: " + (end-start));
*/		
	}

	
	public void OOA(Graph<String,String> graph, SortedSet<String> geneSet, double beta, double delta, String dir, String heatScoreFile) throws IOException{
		PrimitiveMatrix F = createDiffusionMatrixOJA(graph, geneSet, beta);
		PrimitiveMatrix tempExchangedHeatMatrix = createExchangedHeatMatrixOJA(dir, heatScoreFile, F, geneSet);
		RealMatrix exchangedHeatMatrix1 = convertOJAToACM(tempExchangedHeatMatrix);
 		int maxSubnetworkSize1 = obtainMaxSubnetworkSizeForEdgeWeight(exchangedHeatMatrix1, delta);
	}
	
	public void OOO (Graph<String,String> graph, SortedSet<String> geneSet, double beta, double delta, String dir, String heatScoreFile) throws IOException{
		PrimitiveMatrix F2 = createDiffusionMatrixOJA(graph, geneSet, beta);
		PrimitiveMatrix tempExchangedHeatMatrix2 = createExchangedHeatMatrixOJA(dir, heatScoreFile, F2, geneSet);
		int maxSubnetworkSize2 = obtainMaxSubnetworkSizeForEdgeWeight(tempExchangedHeatMatrix2, delta);
	}

	public void OAA (Graph<String,String> graph, SortedSet<String> geneSet, double beta, double delta, String dir, String heatScoreFile) throws IOException{
		PrimitiveMatrix F_2 = createDiffusionMatrixOJA(graph, geneSet, beta);
		RealMatrix F_1 = convertOJAToACM(F_2);
		RealMatrix exchangedHeatMatrix_1 = createExchangedHeatMatrix(dir, heatScoreFile, F_1, geneSet);
		int maxSubnetworkSize = obtainMaxSubnetworkSizeForEdgeWeight(exchangedHeatMatrix_1, delta);
	}
	
	@Test
	public void testACM() throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils(); 

		Path currentPath = Paths.get("");
		String directory = currentPath.toAbsolutePath().toString();
		String fiFile = "FIsInGene_temp.txt";
//		String fiFile = "FIsInGene_031516_with_annotations.txt";
		String heatScoreFile = "heatScore.txt";
		
		//Create the largest component using the whole ReactomeFI network graph
		Graph<String, String> allGenesGraph = gu.createReactomeFIGraph(directory, fiFile);
		Graph<String, String> largestComponent = gu.createLargestComponentGraph(allGenesGraph);
		System.out.println("size: " + largestComponent.getVertices().size());

		//Set of genes in the largest component, ordering determines matrix content ordering
		SortedSet<String> geneSet = gu.getGeneGraphSet(largestComponent);
		
		//HotNet2 Algorithm with Apache Commons Math for all steps
		double beta = 0.5;		//get beta from user
		double delta = 0.05;	//get delta from user
		RealMatrix F = createDiffusionMatrix(largestComponent, geneSet, beta);
//		createTempHeatScoreFile(directory, fiFile, heatScoreFile); //REMOVE THIS and use real heat scores, this is just for testing
		HashMap<String, Double> heatScoreMap = getHeatScoreMap(directory, heatScoreFile, geneSet);
		RealMatrix E = createExchangedHeatMatrix(F, geneSet, heatScoreMap);
		RealMatrix H = identifyHotSubnetworks(E, delta);
		obtainSubnetworkSet(H, geneSet);
//		String diffusionTempFile = "diffusion.txt";
//		saveMatrix(directory+"/output/", diffusionTempFile, F);		
	}

	@Test
	public void testOJA() throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
		String fiFile = "FIsInGene_temp.txt";
//		String fiFile = "FIsInGene_031516_with_annotations.txt";
		String heatScoreFile = "heatScore.txt";
		
		//Create the largest component using the whole ReactomeFI network graph
		Graph<String, String> allGenesGraph = gu.createReactomeFIGraph(directory, fiFile);
		Graph<String, String> largestComponent = gu.createLargestComponentGraph(allGenesGraph);
		System.out.println("size: " + largestComponent.getVertices().size());

		//Set of genes in the largest component, ordering determines matrix content ordering
		SortedSet<String> geneSet = gu.getGeneGraphSet(largestComponent);
		
		//HotNet2 Algorithm with ojAlgo diffusion
		double beta = 0.5;		//get beta from user
		double delta = 0.05;	//get delta from user
		PrimitiveMatrix tempF = createDiffusionMatrixOJA(largestComponent, geneSet, beta);
		RealMatrix F = convertOJAToACM(tempF); 
		HashMap<String, Double> heatScoreMap = getHeatScoreMap(directory, heatScoreFile, geneSet);
		RealMatrix E = createExchangedHeatMatrix(F, geneSet, heatScoreMap);
		RealMatrix H = identifyHotSubnetworks(E, delta);		
		obtainSubnetworkSet(H, geneSet);
	}

	@Test
	public void testOjaToACM() throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
		String fiFile = "FIsInGene_temp.txt";
//		String fiFile = "FIsInGene_031516_with_annotations.txt";
		String heatScoreFile = "heatScore.txt";
		
		//Create the largest component using the whole ReactomeFI network graph
		Graph<String, String> allGenesGraph = gu.createReactomeFIGraph(directory, fiFile);
		Graph<String, String> largestComponent = gu.createLargestComponentGraph(allGenesGraph);
		System.out.println("size: " + largestComponent.getVertices().size());

		//Set of genes in the largest component, ordering determines matrix content ordering
		SortedSet<String> geneSet = gu.getGeneGraphSet(largestComponent);
				
		//HotNet2 Algorithm with ojAlgo diffusion
		double beta = 0.5;	//get beta from user

		//Converts OJA to ACM
		PrimitiveMatrix ojaM = createDiffusionMatrixOJA(largestComponent, geneSet, beta);
		RealMatrix acmM = convertOJAToACM(ojaM);
		System.out.println("ojaM: ------------\n" + ojaM);
		System.out.println("acmM: ------------\n" + acmM);		
	}

	/**
	 * HotNet2 Algorithm
	 * @param directory - Directory location of network interaction file and heat score file.
	 * @param interactionFile - Network interaction file. 
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space.
	 * @param beta - Fraction of own heat each gene retains.  
	 * @param delta - Minimum edge weight threshold; values below edgeWeight are set to 0. 
	 * @throws IOException
	 */
	public void testHotNet2Algorithm(String directory, String interactionFile, String heatScoreFile, double beta, double delta) throws IOException{
		GraphUtils gu = new GraphUtils();
		//Create the largest component using the whole ReactomeFI network graph
		Graph<String, String> allGenesGraph = gu.createReactomeFIGraph(directory, interactionFile);
		Graph<String, String> largestComponent = gu.createLargestComponentGraph(allGenesGraph);

		//Set of genes in the largest component, ordering determines matrix content ordering
		SortedSet<String> geneSet = gu.getGeneGraphSet(largestComponent);
		//HotNet2 Algorithm
		hotnet2Algorithm(directory, heatScoreFile, largestComponent, geneSet, beta, delta);		
	}
	
	/**
	 * Performs the HotNet2 Algorithm
	 * @param directory - Directory where heatScoreFile is located.
	 * @param heatScoreFile - File name containing gene and heat score separated by a space.
	 * @param largestComponent - Graph containing only the largest component in the network.
	 * @param geneSet - Set of genes in the graph used to determine matrix order. 
	 * @param beta - Fraction of own heat each gene retains. 
	 * @param delta - Minimum edge weight threshold; values below edgeWeight are set to 0.
	 * @throws IOException
	 */
	private void hotnet2Algorithm(String directory, String heatScoreFile, Graph<String,String> largestComponent, SortedSet<String> geneSet, double beta, double delta) throws IOException{
		PrimitiveMatrix tempF = createDiffusionMatrixOJA(largestComponent, geneSet, beta);
		RealMatrix F = convertOJAToACM(tempF); 
		HashMap<String, Double> heatScoreMap = getHeatScoreMap(directory, heatScoreFile, geneSet);
		RealMatrix E = createExchangedHeatMatrix(F, geneSet, heatScoreMap);
		RealMatrix H = identifyHotSubnetworks(E, delta);		
		obtainSubnetworkSet(H, geneSet);
	}
	
	/**
	 * Creates an exchanged heat matrix using the provided heat score and beta parameters.
	 * @param directory - Directory where the heatScoreFile is located.
	 * @param heatScoreFile - File name containing gene and heat score separated by a space.
	 * @param graph - Graph containing only the largest component in the network.
	 * @param beta - Fraction of own heat each gene retains.
	 * @return an exchanged heat matrix
	 * @throws IOException
	 */
	public RealMatrix hotNet2ExchangeHeatMatrixFromGraph(String directory, String heatScoreFile, Graph<String,String> graph, double beta) throws IOException{
		GraphUtils gu = new GraphUtils();
		SortedSet<String> geneSet = gu.getGeneGraphSet(graph);
		PrimitiveMatrix F = createDiffusionMatrixOJA(graph, geneSet, beta);
		PrimitiveMatrix tempExchangedHeatMatrix = createExchangedHeatMatrixOJA(directory, heatScoreFile, F, geneSet);
		RealMatrix exchangedHeatMatrix = convertOJAToACM(tempExchangedHeatMatrix);
		return exchangedHeatMatrix;
	}
	
	/**
	 * Creates Diffusion Matrix F=beta(IdentityMatrix-(1-beta)NormalizedAdjacencyMatrix)^(-1).
	 * @param graph - Graph used to create matrix.
	 * @param geneSet - Set of genes in the graph used to determine matrix order.
	 * @param beta - Fraction of own heat each gene retains.
	 * @return a diffusion matrix for HotNet2.
	 */
	public PrimitiveMatrix createDiffusionMatrixOJA(Graph<String, String> graph, SortedSet<String> geneSet, double beta){
		PrimitiveMatrix identityMatrix = PrimitiveMatrix.FACTORY.makeEye(geneSet.size(),geneSet.size());
		double[][] createdAdjMatrix = createNormAdjMatrix(graph, geneSet);
        PrimitiveMatrix adjMatrix = PrimitiveMatrix.FACTORY.rows(createdAdjMatrix);
		PrimitiveMatrix m = adjMatrix.multiply(1-beta);
		m = identityMatrix.subtract(m);
		m = m.invert();
		m = m.multiply(beta);
		return m; 
	}
	
	/**
	 * Creates a normalized adjacency matrix.
	 * @param graph - Graph used to generate adjacency matrix.
	 * @param geneSet - Set of genes used to determine the order and degree of adjacency matrix elements.
	 * @return a normalized adjacency matrix.
	 */
	private double[][] createNormAdjMatrix(Graph<String, String> graph, Set<String> geneSet) {
		GraphUtils gu = new GraphUtils();
		List<String> geneList = new ArrayList<String>(geneSet);
		List<Integer> degreeList = gu.getGeneDegreeList(graph, geneSet);
		int dim = geneList.size();
		double[][] m = new double [dim][dim];
		for (int i=0; i<dim; i++){
			for(int j=0; j<dim; j++){
				String elm1 = geneList.get(i);
				String elm2 = geneList.get(j);				
				if (graph.containsEdge(elm1 + "\t" + elm2) || graph.containsEdge(elm2 + "\t" + elm1)){ 	//CHECK this
					double value = 1.0/degreeList.get(j);
					m[i][j] = value;
				}
			}
		}
		return m; 
	}
	
	/**
	 * Converts ojAlgo PrimitiveMatrix into Apache Commons Math RealMatrix
	 * @param ojaMatrix - PrimitiveMatrix for conversion. 
	 * @return RealMatrix matrix from PrimitiveMatrix.
	 */
	public RealMatrix convertOJAToACM(PrimitiveMatrix ojaMatrix){
		int dim = (int) ojaMatrix.countRows();
		RealMatrix m = MatrixUtils.createRealMatrix(dim, dim);
		for (int i=0; i<dim; i++){
			for(int j=0; j<dim; j++){
				double value = ojaMatrix.get(i, j);
				m.setEntry(i, j, value);
			}
		}
		return m; 
	}	
	
	/**
	 * Creates Diffusion Matrix, F=beta(IdentityMatrix-(1-beta)NormalizedAdjacencyMatrix)^(-1).
	 * @param graph - Graph used to create matrix.
	 * @param geneSet - Set of genes in the graph used to determine matrix order.
	 * @param beta - Fraction of own heat each gene retains.
	 * @return a diffusion matrix for HotNet2.
	 */
	public RealMatrix createDiffusionMatrix (Graph<String, String> graph, Set<String> geneSet, double beta){
		double[][] createdAdjMatrix = createNormAdjMatrix(graph, geneSet);
		RealMatrix m = MatrixUtils.createRealMatrix(createdAdjMatrix);
		m = m.scalarMultiply(1-beta);
		int dim = geneSet.size();
		RealMatrix identityMatrix = MatrixUtils.createRealIdentityMatrix(dim);
		m = identityMatrix.subtract(m);
		m = MatrixUtils.inverse(m);
		m = m.scalarMultiply(beta);
		return m;
	}
	
	/**
	 * Creates Exchanged Heat Matrix, E = F * diagonal heat Matrix.
	 * @param diffusionMatrix - Diffusion matrix created by createDiffusionMatrix().
	 * @param geneSet - Set of genes in the graph used to determine matrix order.
	 * @return a exchanged heat matrix for HotNet2
	 */
	public RealMatrix createExchangedHeatMatrix(String directory, String heatScoreFile, RealMatrix diffusionMatrix, Set<String> geneSet) throws IOException {
		HashMap<String, Double> heatScoreMap = getHeatScoreMap(directory, heatScoreFile, geneSet);
		return createExchangedHeatMatrix(diffusionMatrix, geneSet, heatScoreMap);
	}
	
	/**
	 * Creates Exchanged Heat Matrix, E = F * diagonal heat Matrix, using provided diffusion matrix and heat scores.
	 * @param diffusionMatrix - Diffusion matrix created by createDiffusionMatrix().
	 * @param geneSet - Set of genes in the graph used to determine matrix order.
	 * @param heatScoreMap - HashMap with gene as key and gene's heat score as value
	 * @return a exchanged heat matrix for HotNet2
	 */	
	private RealMatrix createExchangedHeatMatrix(RealMatrix diffusionMatrix, Set<String> geneSet, HashMap<String, Double> heatScoreMap) {
		int dim = diffusionMatrix.getColumnDimension();
		RealMatrix heatMatrix = MatrixUtils.createRealMatrix(dim, dim);
		List<String> geneList = new ArrayList<String>(geneSet);
		for (int i=0; i<dim; i++){
			String gene = geneList.get(i);
			if (heatScoreMap.get(gene) != null)
				heatMatrix.setEntry(i, i, heatScoreMap.get(gene));
		}
		heatMatrix = diffusionMatrix.multiply(heatMatrix);
		return heatMatrix; 
	}
	
	/**
	 * Creates Exchanged Heat Matrix, E = F * diagonal heat Matrix, using provided diffusion matrix and heat scores.
	 * @param directory - Directory location of heat score file.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space.
	 * @param diffusionMatrix - Diffusion matrix created by createDiffusionMatrix().
	 * @param geneSet - Set of genes in the graph used to determine matrix order.
	 * @return a exchanged heat matrix for HotNet2
	 * @throws IOException
	 */
	public PrimitiveMatrix createExchangedHeatMatrixOJA (String directory, String heatScoreFile, PrimitiveMatrix diffusionMatrix, Set<String> geneSet) throws IOException {
		HashMap<String, Double> heatScoreMap = getHeatScoreMap(directory, heatScoreFile, geneSet);
		return createExchangedHeatMatrixOJA(diffusionMatrix, geneSet, heatScoreMap);
	}

	/**
	 * Creates Exchanged Heat Matrix, E = F * diagonal heat Matrix, using provided diffusion matrix and heat scores.
	 * @param diffusionMatrix - Diffusion matrix created by createDiffusionMatrix().
	 * @param geneSet - Set of genes in the graph used to determine matrix order.
	 * @param heatScoreMap - HashMap with gene as key and gene's heat score as value
	 * @return a exchanged heat matrix for HotNet2
	 */
	public PrimitiveMatrix createExchangedHeatMatrixUseMapOJA(PrimitiveMatrix diffusionMatrix, Set<String> geneSet, HashMap<String,Double> heatScoreMap) throws IOException {
		return createExchangedHeatMatrixOJA(diffusionMatrix, geneSet, heatScoreMap);
	}

	/**
	 * Creates Exchanged Heat Matrix, E = F * diagonal heat Matrix, using provided diffusion matrix and heat scores.
	 * @param diffusionMatrix - Diffusion matrix created by createDiffusionMatrix().
	 * @param geneSet - Set of genes in the graph used to determine matrix order.
	 * @param heatScoreMap - HashMap with gene as key and gene's heat score as value
	 * @return a exchanged heat matrix for HotNet2
	 */
	private PrimitiveMatrix createExchangedHeatMatrixOJA(PrimitiveMatrix diffusionMatrix, Set<String> geneSet, HashMap<String, Double> heatScoreMap) {
		int dim = (int) diffusionMatrix.countColumns();
		double[][] tempHeatMatrix = new double [dim][dim];
		List<String> geneList = new ArrayList<String>(geneSet);
		for (int i=0; i<dim; i++){
			String gene = geneList.get(i);
			if (heatScoreMap.get(gene) != null)
				tempHeatMatrix[i][i] = heatScoreMap.get(gene);
		}		
		PrimitiveMatrix heatMatrix = PrimitiveMatrix.FACTORY.rows(tempHeatMatrix);
		heatMatrix = diffusionMatrix.multiplyRight(heatMatrix);
		
		return heatMatrix; 
	}
	
	/**
	 * Creates Exchanged Heat Matrix, E = F * diagonal heat Matrix, using provided diffusion matrix and heat scores.
	 * @param diffusionMatrix - Diffusion matrix created by createDiffusionMatrix().
	 * @param geneSet - Set of genes in the graph used to determine matrix order.
	 * @param heatScoreMap - HashMap with gene as key and gene's heat score as value
	 * @return a exchanged heat matrix for HotNet2
	 */
	public RealMatrix createExchangedHeatMatrixUseMap(RealMatrix diffusionMatrix, Set<String> geneSet, HashMap<String,Double> heatScoreMap) throws IOException {
		int dim = diffusionMatrix.getRowDimension();
		RealMatrix heatMatrix = MatrixUtils.createRealMatrix(dim, dim);
		List<String> geneList = new ArrayList<String>(geneSet);
		for (int i=0; i<dim; i++){
			String gene = geneList.get(i);
			if (heatScoreMap.get(gene) != null)
				heatMatrix.addToEntry(i, i, heatScoreMap.get(gene));
		}		
		heatMatrix = diffusionMatrix.multiply(heatMatrix); 

		return heatMatrix;
	}
	
	/**
	 * Obtain a set of unique exchanged heat matrix values using RealMatrix. 
	 * @param matrix - Matrix used to obtain element values.
	 * @return a sorted set of element values.
	 */
	public SortedSet<Double> obtainUniqueExchangedHeatMatrixValues (RealMatrix matrix){
		SortedSet<Double> weights = new TreeSet<Double>();
		int dim = matrix.getRowDimension();		
		for (int i=0; i<dim; i++){
			for(int j=0; j<dim; j++)
				weights.add(matrix.getEntry(i, j));
		}
		return weights;
	}

	/**
	 * Obtain a set of unique exchanged heat matrix values using PrimitiveMatrix. 
	 * @param matrix - Matrix used to obtain element values.
	 * @return a sorted set of element values.
	 */
	public SortedSet<Double> obtainUniqueExchangedHeatMatrixValues (PrimitiveMatrix matrix){
		SortedSet<Double> weights = new TreeSet<Double>();
		int dim = (int) matrix.countRows();	
		for (int i=0; i<dim; i++){
			for(int j=0; j<dim; j++)
				weights.add(matrix.get(i, j));
		}
		return weights;
	}
	
	/**
	 * Stores heat scores of genes from a file into a HashMap, only genes in the provided set are included. 
	 * @param directory - Directory of file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space. 
	 * @param geneSet - Set of genes that should be included in result.
	 * @return a HashMap with the gene as key and it's heat score as value.
	 * @throws IOException
	 */
	public HashMap<String, Double> getHeatScoreMap(String directory, String heatScoreFile, Set<String> geneSet) throws IOException {
		HashMap<String, Double> geneScores = new HashMap<String, Double>();
		Path filePath = Paths.get(directory, heatScoreFile);
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(filePath, charset);		
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(" ");
			String gene = row[0];
			if (geneSet.contains(gene)){
				Double score = Double.parseDouble(row[1]);
				geneScores.put(gene, score);
			}
		}
		br.close();
		
		Set<String> geneWithScore = new HashSet<String>(geneScores.keySet());
		for (String gene: geneSet){
			if(!geneWithScore.contains(gene))
				geneScores.put(gene, 0.0);
		}	
		return geneScores;
	}
	
	/**
	 * Gets genes with missing heat scores.
	 * @param directory - Directory of file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space.
	 * @param graph - Graph used to obtain set of genes to find heat scores.
	 * @throws IOException 
	 */
	public void getMissingHeatGenes(String directory, String heatScoreFile, Graph<String,String> graph) throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		String geneMissingHeatFile = "geneMissingHeat.txt";
		SortedSet<String> geneSet = gu.getGeneGraphSet(graph);
		Set<String> geneMissingHeat = getMissingHeatScoreGenes(directory, heatScoreFile, geneSet);
		fu.saveSetToFile(directory+"/output/", geneMissingHeatFile, geneMissingHeat);
	}

	/**
	 * Gets genes with missing heat scores by comparing heat score file genes with genes from a set.
	 * @param directory - Directory of file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space.
	 * @param geneSet - Set of genes that should be included in result.
	 * @return a set of genes and scores ordered by geneSet.
	 * @throws IOException
	 */
	private Set<String> getMissingHeatScoreGenes(String directory, String heatScoreFile, Set<String> geneSet) throws IOException {
		Set<String> heatScoreGenes = new HashSet<String>();
		Path filePath = Paths.get(directory, heatScoreFile);
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(filePath, charset);
		br.readLine();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(" ");
			String gene = row[0];
			heatScoreGenes.add(gene);
		}
		br.close();
		geneSet.removeAll(heatScoreGenes);
		return geneSet;
	}
		
	/**
	 * Creates a file containing random heat scores for each gene in the ReactomeFI network
	 * <p>
	 * <b>Note:</b> this is for testing purposes only.  
	 * @param directory	- Directory where files are located and will be stored.
	 * @param fiFile - Reactome FI network file name.
	 * @param heatScoreFile - File name of genes and their heat score. 
	 * @throws IOException
	 */
	public void createTempHeatScoreFile(String directory, String fiFile, String heatScoreFile) throws IOException{
		FileUtils fu = new FileUtils();
		Path fiFilePath = Paths.get(directory, fiFile);
		Set<String> allGenes = fu.getAllGenesReactome(fiFilePath);
		String filePath = directory + "/" + heatScoreFile;
		PrintWriter pw = new PrintWriter(filePath);
		Random randomGen = new Random();
		for(String g: allGenes){
			Double random = randomGen.nextDouble();
			pw.println(g + " " + random);
		}
		pw.close();
	}
	
	/** 
	 * Identifies Hot Subnetworks by removing values in exchanged heat matrix below the provided threshold. 
	 * <p>
	 * <b>Note:</b> modified Apache Commons Math3's walkInOptimizedOrder() to replace values below the threshold as 0.
	 * @param matrix - Exchanged heat matrix.
	 * @param edgeWeight - Minimum edge weight threshold; values below edgeWeight are set to 0. 
	 * @return a matrix with values below edgeWeight removed
	 */
	public RealMatrix identifyHotSubnetworks(RealMatrix matrix, double edgeWeight){
		final double weight = edgeWeight;
		RealMatrix exchangedHeatMatrix = matrix.copy();
		exchangedHeatMatrix.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
			public double visit(int row, int column, double value){	
				if (value < weight)
					return 0;	
				return value; 
			}
		});	
		return exchangedHeatMatrix;
	}
	
	/**
	 * Identifies size of the largest Hot Subnetworks by removing values in exchanged heat matrix below the provided threshold. 
	 * @param matrix - Exchanged heat matrix.
	 * @param edgeWeight - minimum edge weight threshold; values below edgeWeight are set to 0.
	 * @return the size of the largest subnetwork.
	 */
	public int identifyHotSubnetworks2(RealMatrix matrix, double edgeWeight){
		final double weight = edgeWeight;
		RealMatrix exchangedHeatMatrix = matrix.copy();
		exchangedHeatMatrix.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
			public double visit(int row, int column, double value){	
				if (value < weight)
					return 0;	
				return value; 
			}
		});
		DefaultDirectedGraph<Integer, DefaultEdge> graph = new DefaultDirectedGraph<Integer, DefaultEdge>(DefaultEdge.class);
		int size = exchangedHeatMatrix.getRowDimension();
		//Add vertices to the graph
		for (int i=0; i<size; i++)
			graph.addVertex(i);
		//Add edges greater than and equal the edgeweight to the graph 
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (i != j && exchangedHeatMatrix.getEntry(i, j) >= edgeWeight)
						graph.addEdge(i,j);
			}
		}
		StrongConnectivityInspector sci = new StrongConnectivityInspector(graph);
		List<Set<Integer>> components = sci.stronglyConnectedSets();
		int subnetworkSize = 0;
		for(Set<Integer> c: components){
			if(c.size() > subnetworkSize)
				subnetworkSize = c.size();
		}
		return subnetworkSize;
	}
	
	/**
	 * Obtain max subnetwork size for the provided edge weight.
	 * <p>
	 * <b>Note:</> the graph generated from the matrix is a directed graph and subnetworks obtained are strongly connected. 
	 * @param matrix - Matrix to extract subnetworks from.
	 * @param edgeWeight - Edge weight value is used to remove edges less than this value to create strongly connected components.
	 * @return  the size of the largest component.
	 */
	public int obtainMaxSubnetworkSizeForEdgeWeight (RealMatrix matrix, double edgeWeight){
		DefaultDirectedGraphModified<Integer, DefaultEdge> graph = new DefaultDirectedGraphModified<Integer, DefaultEdge>(DefaultEdge.class);
		int size = matrix.getRowDimension();
		//Add vertices to the graph
		for (int i=0; i<size; i++)
			graph.addVertex(i);
		//Add edges greater than and equal the edgeweight to the graph 
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (i != j && matrix.getEntry(i, j) >= edgeWeight)
						graph.addEdge(i,j);
			}
		}				
		StrongConnectivityInspector sci = new StrongConnectivityInspector(graph);
		List<Set<Integer>> components = sci.stronglyConnectedSets();
		int subnetworkSize = 0;
		for(Set<Integer> c: components){
			if(c.size() > subnetworkSize)
				subnetworkSize = c.size();
		}
		return subnetworkSize;
	}

	/**
	 * Obtain max subnetwork size for the provided edge weight.(Uses JUNG and code written by chschmitz)
	 * <p>
	 * <b>Note:</> the graph generated from the matrix is a directed graph and subnetworks obtained are strongly connected. 
	 * @param matrix - Matrix to extract subnetworks from.
	 * @param edgeWeight - Edge weight value is used to remove edges less than this value to create strongly connected components.
	 * @return  the size of the largest component.
	 */	
	public int obtainMaxSubnetworkSizeForEdgeWeightJUNG (RealMatrix exchangedHeatMatrix, double edgeWeight){
		DirectedGraph<Integer, String> graph = new DirectedSparseGraph<Integer, String>();
		int size = exchangedHeatMatrix.getRowDimension();
		//Add vertices to the graph
		for (int i=0; i<size; i++)
			graph.addVertex(i);
		//Add edges greater than and equal the edgeweight to the graph 
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (i != j && exchangedHeatMatrix.getEntry(i, j) >= edgeWeight)
						graph.addEdge(i + "\t" + j, i,j);
			}
		}
		Collection<Set<Integer>> components = StronglyConnectedComponents.strongComponentsAsSets(graph);	
		int subnetworkSize = 0;
		for(Set<Integer> c: components){
			if(c.size() > subnetworkSize)
				subnetworkSize = c.size();
		}
		return subnetworkSize;
	}
	
	/**
	 * Obtain the max subnetwork size from a matrix using edge weight as a threshold.
	 * @param matrix - Matrix to extract subnetworks from. 
	 * @param edgeWeight - Edge weight value used to remove edges less than this value to create strongly connected components.
	 * @return the maximum subnetwork size.
	 */
	public int obtainMaxSubnetworkSizeForEdgeWeight (PrimitiveMatrix matrix, double edgeWeight){
		DefaultDirectedGraph<Integer, DefaultEdge> graph = new DefaultDirectedGraph<Integer, DefaultEdge>(DefaultEdge.class);
		int size = (int) matrix.countRows();
		//Add vertices to the graph
		for (int i=0; i<size; i++)
			graph.addVertex(i);
		//Add edges greater than and equal the edgeweight to the graph 
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (i != j && matrix.get(i, j) >= edgeWeight)
						graph.addEdge(i,j);
			}
		}				
		StrongConnectivityInspector sci = new StrongConnectivityInspector(graph);
		List<Set<Integer>> components = sci.stronglyConnectedSets();
		int subnetworkSize = 0;
		for(Set<Integer> c: components){
			if(c.size() > subnetworkSize)
				subnetworkSize = c.size();
		}
		return subnetworkSize;
	}

	/**
	 * Obtain a directed graph from a matrix using edge weight as a threshold.
	 * @param matrix - Matrix to extract subnetworks from. 
	 * @param edgeWeight - Edge weight value used to remove edges less than this value to create strongly connected components.
	 * @return a directed graph. 
	 */
	public DefaultDirectedGraph<String, DefaultEdge> obtainGraphForEdgeWeight (RealMatrix matrix, double edgeWeight){
		DefaultDirectedGraph<String, DefaultEdge> graph = new DefaultDirectedGraph<String, DefaultEdge>(DefaultEdge.class);
		int size = matrix.getRowDimension();
		//Add vertices to the graph
		for (int i=0; i<size; i++)
			graph.addVertex(Integer.toString(i));
		//Add edges greater than and equal the edgeweight to the graph 
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (i != j && matrix.getEntry(i, j) >= edgeWeight)
						graph.addEdge(Integer.toString(i),Integer.toString(j));
			}
		}
		return graph;
	}
	
	/**
	 * Obtain largest subnetwork size from the provided matrix.  
	 * @param matrix - Matrix used to determine subnetwork sizes.
	 * @param geneSet - Set of genes used to determine matrix order. 
	 * @return the value of the largest component size.
	 */
	public int obtainLargestSubnetworkSize(RealMatrix matrix, SortedSet<String> geneSet){
		GraphUtils gu = new GraphUtils();
		RealMatrix tempMatrix = matrix.copy(); 
		Graph<String, String> graph = gu.covertMatrixToDirectedGraph(tempMatrix, geneSet);
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		int largestSize = 0;
		for (Set<String> c: components){
			if(c.size() > largestSize)
				largestSize = c.size();
		}
		return largestSize;
	}
	
	/**
	 * Obtains a set of subnetwork genes.
	 * @param matrix - Matrix to extract subnetworks from.
	 * @param geneSet - Set of genes used to determine matrix order.
	 * @return a set containing sets of genes that form subnetworks.
	 */
	public Set<Set<String>> obtainSubnetworkSet(RealMatrix matrix, SortedSet<String> geneSet){
		GraphUtils gu = new GraphUtils();
		Graph<String, String> graph = gu.covertMatrixToDirectedGraph(matrix, geneSet);
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		
		System.out.println("There were " + components.size() + " subnetwork(s) identified");
	
		List<Set<String>> componentList = new ArrayList<Set<String>>(components);	
		for (int i=0; i<components.size(); i++){
			List<String> component = new ArrayList(componentList.get(i));
			String gene = component.get(0);
			
			graph.getOutEdges(gene);//Edges for a specific component
			//********** save each subnetwork into a edge list file?
		}		
		return components;
	}

	/**
	 * Saves specified RealMatrix into a tab delimited file.
	 * @param directory - Directory to store file.
	 * @param fileName - Name of file to store matrix contents.
	 * @param matrix a RealMatrix to be stored.
	 * @throws IOException
	 */
	public void saveMatrix(String directory, String fileName, RealMatrix matrix) throws IOException{
		File file = new File(directory+"/output/"+ fileName);
		file.getParentFile().mkdirs();
		if (!file.exists())
			file.createNewFile();
		String filePath = directory + "/output/" + fileName;
		PrintWriter pw = new PrintWriter(filePath);
		int matrixDim = matrix.getRowDimension();
		for (int i=0; i<matrixDim; i++){
			double[] rowValues = matrix.getRow(i);
			for (double r: rowValues)
				pw.print(r + "\t");
			pw.println();
		}
		pw.close();
	}
		 
	/**
	 * Saves specified PrimitiveMatrix into a tab delimited file.
	 * @param directory - Directory to store file.
	 * @param fileName- File name to store matrix contents.
	 * @param matrix - PrimitiveMatrix to be stored.
	 * @throws IOException
	 */
	public void saveMatrix(String directory, String fileName, PrimitiveMatrix matrix) throws IOException{
		File file = new File(directory+"/output/"+ fileName);
		file.getParentFile().mkdirs();
		if (!file.exists())
			file.createNewFile();
		String filePath = directory + "/output/" + fileName;
		PrintWriter pw = new PrintWriter(filePath);
		int matrixDim = (int) matrix.countRows();
		for (int i=0; i<matrixDim; i++){
			for (int j=0; j<matrixDim; j++){
				pw.print(matrix.get(i,j) + "\t");
			}
			pw.println();
		}
		pw.close();
	}
	
}