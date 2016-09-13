package hotnet2;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Test;
import org.ojalgo.matrix.PrimitiveMatrix;

import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.Graph;

public class HotNet2Matrices {
	
	public HotNet2Matrices(){
	
	}
	
	@Test
	public void TestACM() throws IOException{
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
		String fiFile = "FIsInGene_temp.txt";
//		String fiFile = "FIsInGene_031516_with_annotations.txt";
		String heatScoreFile = "heatScore.txt";
		
		//Create the largest component using the whole ReactomeFI network graph
		GraphUtils gu = new GraphUtils();
		Graph<String, String> allGenesGraph = gu.createReactomeFIGraphWrapper(directory, fiFile);
		Graph<String, String> largestComponent = gu.createLargestComponentGraphWrapper(allGenesGraph);
		System.out.println("size: " + largestComponent.getVertices().size());

		//Set of genes in the largest component, ordering determines matrix content ordering
		Set<String> geneSet = gu.getGeneGraphSetWrapper(largestComponent);
		
		//HotNet2 Algorithm with Apache Commons Math for all steps
		BigDecimal tempBeta = new BigDecimal("0.5");	//get beta from user
		double beta = tempBeta.doubleValue();
		BigDecimal tempDelta = new BigDecimal("0.05");	//get delta from user
		double delta = tempDelta.doubleValue();
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
	public void TestOJA() throws IOException{
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
		String fiFile = "FIsInGene_temp.txt";
//		String fiFile = "FIsInGene_031516_with_annotations.txt";
		String heatScoreFile = "heatScore.txt";
		
		//Create the largest component using the whole ReactomeFI network graph
		GraphUtils gu = new GraphUtils();
		Graph<String, String> allGenesGraph = gu.createReactomeFIGraphWrapper(directory, fiFile);
		Graph<String, String> largestComponent = gu.createLargestComponentGraphWrapper(allGenesGraph);
		System.out.println("size: " + largestComponent.getVertices().size());

		//Set of genes in the largest component, ordering determines matrix content ordering
		Set<String> geneSet = gu.getGeneGraphSetWrapper(largestComponent);
		
		//HotNet2 Algorithm with ojAlgo diffusion
		BigDecimal tempBeta = new BigDecimal("0.5");	//get beta from user
		double beta = tempBeta.doubleValue();
		BigDecimal tempDelta = new BigDecimal("0.05");	//get delta from user
		double delta = tempDelta.doubleValue();
		PrimitiveMatrix tempF = createDiffusionMatrixOJA(largestComponent, geneSet, beta);
		RealMatrix F = convertOJAToACM(tempF); 
		HashMap<String, Double> heatScoreMap = getHeatScoreMap(directory, heatScoreFile, geneSet);
		RealMatrix E = createExchangedHeatMatrix(F, geneSet, heatScoreMap);
		RealMatrix H = identifyHotSubnetworks(E, delta);		
		obtainSubnetworkSet(H, geneSet);			
	}

	@Test
	public void TestOjaToACM() throws IOException{
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
		String fiFile = "FIsInGene_temp.txt";
//		String fiFile = "FIsInGene_031516_with_annotations.txt";
		String heatScoreFile = "heatScore.txt";
		
		//Create the largest component using the whole ReactomeFI network graph
		GraphUtils gu = new GraphUtils();
		Graph<String, String> allGenesGraph = gu.createReactomeFIGraphWrapper(directory, fiFile);
		Graph<String, String> largestComponent = gu.createLargestComponentGraphWrapper(allGenesGraph);
		System.out.println("size: " + largestComponent.getVertices().size());

		//Set of genes in the largest component, ordering determines matrix content ordering
		Set<String> geneSet = gu.getGeneGraphSetWrapper(largestComponent);
				
		//HotNet2 Algorithm with ojAlgo diffusion
		BigDecimal tempBeta = new BigDecimal("0.5");	//get beta from user
		double beta = tempBeta.doubleValue();

		//Converts OJA to ACM
		PrimitiveMatrix ojaM = createDiffusionMatrixOJA(largestComponent, geneSet, beta);
		RealMatrix acmM = convertOJAToACM(ojaM);
		System.out.println("ojaM: ------------\n" + ojaM);
		System.out.println("acmM: ------------\n" + acmM);		
	}

	
	/**
	 * Creates Diffusion Matrix F=beta(IdentityMatrix-(1-beta)NormalizedAdjacencyMatrix)^(-1).
	 * @param graph - Graph used to create matrix.
	 * @param geneSet - Set of genes in the graph used to determine matrix order.
	 * @param beta - Fraction of own heat each gene retains.
	 * @return a diffusion matrix for HotNet2.
	 */
	public PrimitiveMatrix createDiffusionMatrixOJAWrapper(Graph<String, String> graph, Set<String> geneSet, double beta){
		return createDiffusionMatrixOJA(graph, geneSet, beta);
	}
	
	/**
	 * Creates Diffusion Matrix F=beta(IdentityMatrix-(1-beta)NormalizedAdjacencyMatrix)^(-1).
	 * @param graph - Graph used to create matrix.
	 * @param geneSet - Set of genes in the graph used to determine matrix order.
	 * @param beta - Fraction of own heat each gene retains.
	 * @return a diffusion matrix for HotNet2.
	 */
	private PrimitiveMatrix createDiffusionMatrixOJA(Graph<String, String> graph, Set<String> geneSet, double beta){
		PrimitiveMatrix identityMatrix = PrimitiveMatrix.FACTORY.makeEye(geneSet.size(),geneSet.size());
		double[][] createdAdjMatrix = createNormAdjMatrixOJA(graph, geneSet);
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
	private double[][] createNormAdjMatrixOJA(Graph<String, String> graph, Set<String> geneSet) {
		GraphUtils gu = new GraphUtils(); 
		List<String> geneList = new ArrayList<String>(geneSet);
		List<Integer> degreeList = gu.getGeneDegreeListWrapper(graph, geneSet);		
		int dim = geneList.size();
		double[][] m = new double [dim][dim];
		for (int i=0; i<dim; i++){
			for(int j=0; j<dim; j++){
				String elm1 = geneList.get(i);
				String elm2 = geneList.get(j);				
				if (graph.containsEdge(elm1 + "\t" + elm2) || graph.containsEdge(elm2 + "\t" + elm1)){ 	//CHECK this
					double value = 1.0/degreeList.get(i);
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
	public RealMatrix convertOJAToACMWrapper(PrimitiveMatrix ojaMatrix){
		return convertOJAToACM(ojaMatrix);
	}
	
	/**
	 * Converts ojAlgo PrimitiveMatrix into Apache Commons Math RealMatrix
	 * @param ojaMatrix - PrimitiveMatrix for conversion. 
	 * @return RealMatrix matrix from PrimitiveMatrix.
	 */
	private RealMatrix convertOJAToACM(PrimitiveMatrix ojaMatrix){
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
	public RealMatrix createDiffusionMatrixWrapper(Graph<String, String> graph, Set<String> geneSet, double beta){
		return createDiffusionMatrix(graph, geneSet, beta);
	}	
	
	/**
	 * Creates Diffusion Matrix, F=beta(IdentityMatrix-(1-beta)NormalizedAdjacencyMatrix)^(-1).
	 * @param graph - Graph used to create matrix.
	 * @param geneSet - Set of genes in the graph used to determine matrix order.
	 * @param beta - Fraction of own heat each gene retains.
	 * @return a diffusion matrix for HotNet2.
	 */
	private RealMatrix createDiffusionMatrix (Graph<String, String> graph, Set<String> geneSet, double beta){
		RealMatrix m = createNormAdjMatrix(graph, geneSet);	
		m = m.scalarMultiply(1-beta);
		int dim = geneSet.size();
		RealMatrix identityMatrix = MatrixUtils.createRealIdentityMatrix(dim);
		m = identityMatrix.subtract(m);
		m = MatrixUtils.inverse(m);
		m = m.scalarMultiply(beta);
		return m;
	}
	
	/**
	 * Creates a normalized adjacency matrix.
	 * @param graph - Graph used to generate adjacency matrix.
	 * @param geneSet - Set of genes used to determine the order and degree of adjacency matrix elements.
	 * @return a normalized adjacency matrix.
	 */
	private RealMatrix createNormAdjMatrix(Graph<String, String> graph, Set<String> geneSet) {
		GraphUtils gu = new GraphUtils();
		
		List<String> geneList = new ArrayList<String>(geneSet);
		List<Integer> degreeList = gu.getGeneDegreeListWrapper(graph, geneSet);		
		int dim = geneList.size();
		RealMatrix m = MatrixUtils.createRealMatrix(dim, dim);
		for (int i=0; i<dim; i++){
			for(int j=0; j<dim; j++){
				String elm1 = geneList.get(i);
				String elm2 = geneList.get(j);				
				if (graph.containsEdge(elm1 + "\t" + elm2) || graph.containsEdge(elm2 + "\t" + elm1)){ 	//CHECK this
					double value = 1.0/degreeList.get(i);
					m.setEntry(i,j,value);
				}
			}
		}
		return m; 
	}

	/**
	 * Creates Exchanged Heat Matrix, E = F * diagonal heat Matrix.
	 * @param diffusionMatrix - Diffusion matrix created by createDiffusionMatrix().
	 * @param geneSet - Set of genes in the graph used to determine matrix order.
	 * @return a exchanged heat matrix for HotNet2
	 */
	public RealMatrix createExchangedHeatMatrixWrapper(String directory, String heatScoreFile, RealMatrix diffusionMatrix, Set<String> geneSet) throws IOException {
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
	 * Stores heat scores of genes from a file into a HashMap, only genes in the provided set are included. 
	 * @param directory - Directory of file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space. 
	 * @param geneSet - Set of genes that should be included in result.
	 * @return a HashMap with the gene as key and it's heat score as value.
	 * @throws IOException
	 */
	public HashMap<String, Double> getHeatScoreMapWrapper(String directory, String heatScoreFile, Set<String> geneSet) throws IOException {
		return getHeatScoreMap(directory, heatScoreFile, geneSet);
	}
	
	/**
	 * Stores heat scores of genes from a file into a HashMap, only genes in the provided set are included. 
	 * @param directory - Directory of file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space. 
	 * @param geneSet - Set of genes that should be included in result.
	 * @return a HashMap with the gene as key and it's heat score as value.
	 * @throws IOException
	 */
	private HashMap<String, Double> getHeatScoreMap(String directory, String heatScoreFile, Set<String> geneSet) throws IOException {
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
		return geneScores;
	}
	
	/**
	 * Gets genes with missing heat scores.
	 * @param directory - Directory of file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space.
	 * @param graph - Graph used to obtain set of genes to find heat scores.
	 * @throws IOException 
	 */
	public void getMissingHeatGenesWrapper(String directory, String heatScoreFile, Graph<String,String> graph) throws IOException{
		getMissingHeatGenes(directory, heatScoreFile, graph);
	}
	
	/**
	 * Gets genes with missing heat scores.
	 * @param directory - Directory of file location.
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space.
	 * @param graph - Graph used to obtain set of genes to find heat scores.
	 * @throws IOException 
	 */
	private void getMissingHeatGenes(String directory, String heatScoreFile, Graph<String,String> graph) throws IOException{
		GraphUtils gu = new GraphUtils();
		FileUtils fu = new FileUtils();
		String geneMissingHeatFile = "geneMissingHeat.txt";
		Set<String> geneSet = gu.getGeneGraphSetWrapper(graph);
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
	public void createTempHeatScoreFileWrapper(String directory, String fiFile, String heatScoreFile) throws IOException{
		createTempHeatScoreFile(directory, fiFile, heatScoreFile);
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
	private void createTempHeatScoreFile(String directory, String fiFile, String heatScoreFile) throws IOException{
		FileUtils fu = new FileUtils();
		Path fiFilePath = Paths.get(directory, fiFile);
		Set<String> allGenes = fu.getAllGenesReactomeWrapper(fiFilePath);
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
	 * @param exchangedHeatMatrix - Exchanged heat matrix.
	 * @param edgeWeight - Minimum edge weight threshold; values below edgeWeight are set to 0. 
	 * @return a matrix with values below edgeWeight removed
	 */
	public RealMatrix identifyHotSubnetworksWrapper(RealMatrix exchangedHeatMatrix, double edgeWeight){
		return identifyHotSubnetworks(exchangedHeatMatrix, edgeWeight);
	}
	
	/** 
	 * Identifies Hot Subnetworks by removing values in exchanged heat matrix below the provided threshold. 
	 * <p>
	 * <b>Note:</b> modified Apache Commons Math3's walkInOptimizedOrder() to replace values below the threshold as 0.
	 * @param exchangedHeatMatrix - Exchanged heat matrix.
	 * @param edgeWeight - Minimum edge weight threshold; values below edgeWeight are set to 0. 
	 * @return a matrix with values below edgeWeight removed
	 */
	private RealMatrix identifyHotSubnetworks(RealMatrix exchangedHeatMatrix, double edgeWeight){
		final double weight= edgeWeight;
		exchangedHeatMatrix.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
			public double visit(int row, int column, double value){	
				if (value <= weight)
					return 0;	
				return value; 
			}
		});	
		return exchangedHeatMatrix;
	}

	/**
	 * Obtains a set of subnetwork genes.
	 * @param matrix - Matrix to extract subnetworks from.
	 * @param geneSet - Set of genes used to determine matrix order.
	 * @return a set containing sets of genes that form subnetworks.
	 */
	private Set<Set<String>> obtainSubnetworkSet (RealMatrix matrix, Set<String> geneSet){
		GraphUtils gu = new GraphUtils(); 
		Graph<String, String> graph = gu.covertMatrixToGraphWrapper(matrix, geneSet);

		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graph);
		System.out.println("There was " + components.size() + " subnetwork(s) identified");
	
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
	private void saveMatrixWrapper(String directory, String fileName, RealMatrix matrix) throws IOException{
		saveMatrix(directory, fileName, matrix);
	}
	
	/**
	 * Saves specified RealMatrix into a tab delimited file.
	 * @param directory - Directory to store file.
	 * @param fileName - Name of file to store matrix contents.
	 * @param matrix a RealMatrix to be stored.
	 * @throws IOException
	 */
	private void saveMatrix(String directory, String fileName, RealMatrix matrix) throws IOException{
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
	public void saveMatrixWrapper(String directory, String fileName, PrimitiveMatrix matrix) throws IOException {
		saveMatrix(directory, fileName, matrix);
	}
	 
	/**
	 * Saves specified PrimitiveMatrix into a tab delimited file.
	 * @param directory - Directory to store file.
	 * @param fileName- File name to store matrix contents.
	 * @param matrix - PrimitiveMatrix to be stored.
	 * @throws IOException
	 */
	private void saveMatrix(String directory, String fileName, PrimitiveMatrix matrix) throws IOException{
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