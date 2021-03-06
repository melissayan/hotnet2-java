package edu.ohsu.hotnet2;

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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Test;
import org.ojalgo.matrix.PrimitiveMatrix;

import edu.uci.ics.jung.algorithms.importance.BetweennessCentrality;
import edu.uci.ics.jung.graph.Graph;

public class BetaSelection{
	private Path currentPath;
	private String directory;
	private String betweennessScoreFile;
	private String influenceFile;
		
	public BetaSelection(){
		this.currentPath = Paths.get("");
		this.directory = currentPath.toAbsolutePath().toString();
		this.betweennessScoreFile = "betweennessScore.txt";
		this.influenceFile = "influence.txt";
	}
	
	@Test
	public void testBetaSelection() throws IOException{
		String fiFile = "FIsInGene_031516_with_annotations.txt";
		Path fiFilePath = Paths.get(directory, fiFile);
		Boolean onlyTP53 = false; 
		selectBeta(directory, fiFile, betweennessScoreFile, influenceFile, onlyTP53);
	}
	
	/**
	 * Selects the beta parameter to assign an amount of heat retained by each gene for creating the diffusion matrix.
	 * <p>
	 * <b>Note:</b> results are saved in a textfile for processing in excel.
	 * @param directory - Directory of Reactome FI network file and place to save files.
	 * @param fiFile - Reactome FI network file.
	 * @param betweennessScoreFile - File name for gene betweenness centrality scores. 
	 * @param influenceFile - File name for saving influence gene counts.
	 * @throws IOException
	 */
	public void selectBeta(String directory, String fiFile, String betweennessScoreFile, String influenceFile, boolean onlyTP53) throws IOException{
		GraphUtils gu = new GraphUtils();
		//Create the largest component using the whole ReactomeFI network graph
		Graph<String, String> allGenesGraph = gu.createReactomeFIGraph(directory, fiFile);
		Graph<String, String> largestComponent = gu.createLargestComponentGraph(allGenesGraph);
		if (onlyTP53 == true){
			//Get TP53 as a source protein
			List<String> sourceProteins =  new ArrayList<String>();
			sourceProteins.add("TP53");
			String newDirectory = directory + "/output/";
			generateRangeBeta(newDirectory, influenceFile, largestComponent, sourceProteins);
		}
		else {
			//Calculate and save betweenness centrality for all genes in largest component
			saveBetweennessCentrality(directory, betweennessScoreFile, largestComponent);

			//Get 5 source proteins from betweennness centrality scores
			String newDirectory = directory + "/output/"; 
			List<String> sourceProteins = getSourceProteins(newDirectory, betweennessScoreFile);
			generateRangeBeta(newDirectory, influenceFile, largestComponent, sourceProteins);
		}
	}
	
	/**
	 * Generates 20 diffusion matrices with different beta parameters ranging from 0.05 to 1.00.
	 * @param directory - Directory of Reactome FI network file and place to save files.
	 * @param influenceFile - File name for saving influence gene counts.
	 * @param largestComponent - Largest component in graph. 
	 * @param influenceFile - File name for saving influence gene counts.
	 * @throws IOException
	 */
	private void generateRangeBeta(String directory, String influenceFile, Graph<String, String> largestComponent, List<String> sourceProteins) throws IOException{
		GraphUtils gu = new GraphUtils();
		HotNet2Matrix hn2m = new HotNet2Matrix();
		//Generate 20 diffusion matrices with different beta ranging from: 0.05, 0.10, 0.15,..., 1.00
		SortedSet<String> geneSet = gu.getGeneGraphSet(largestComponent);
		for(int i=1; i<21; i++){
			BigDecimal tempBeta = new BigDecimal("0.05");
			tempBeta = tempBeta.multiply(new BigDecimal(i));
			double beta = tempBeta.doubleValue();
			System.out.println("----beta: "+ beta);
			long start1 = System.currentTimeMillis();
			PrimitiveMatrix diffusionMatrix = hn2m.createDiffusionMatrixOJA(largestComponent, geneSet, beta);
			long end1 = System.currentTimeMillis();
			System.out.println("\tDiffusion Matrix Time Taken: " + ((end1 - start1) / 1000) + " seconds");
//			String newDirectory = directory +"/output/inflectionPoint"; 
//			String newFileName = tempBeta + "_diffusion.txt";
//			saveMatrix(newDirectory, newFileName, diffusionMatrix);

			saveInfluenceGeneCountRange(directory, influenceFile, largestComponent, sourceProteins, diffusionMatrix, geneSet, tempBeta);			
		}		
	}
	
	
	/**
	 * Selects the beta parameter for the iRefIndex network.
	 * <p>
	 * <b>Note:</b> to compare implementation matches HotNet2 Supplementary Figure 24 for iRefIndex network with beta=0.45. 
	 * @param directory - Directory of iRefIndex network file and place to save files. 
	 * @param file - iRefIndex network file.
	 * @param betweennessScoreFile - File name for gene betweenness centrality scores. 
	 * @param influenceFile - File name for saving influence gene counts.
	 * @throws IOException
	 */
	public void selectBetaForIrefindex(String directory, String file, String betweennessScoreFile, String influenceFile) throws IOException{
		GraphUtils gu = new GraphUtils();
		HotNet2Matrix hn2m = new HotNet2Matrix(); 

		//Create the largest component using the iRefIndex network graph
		Graph<String, String> allGenesGraph = gu.createReactomeFIGraph(directory, file);
		Graph<String, String> largestComponent = gu.createLargestComponentGraph(allGenesGraph);
		
		//Get TP53 source protein from iref_edge_list
		List<String> sourceProteins =  new ArrayList<String>();
		sourceProteins.add("10922");
//		sourceProteins.add("6911"); //source protein for iref_edge_list_temp
  
		//Generate diffusion matrix with beta: 0.45
		SortedSet<String> geneSet = gu.getGeneGraphSet(largestComponent);
		for(int i=1; i<2; i++){
			BigDecimal tempBeta = new BigDecimal("0.45");
			tempBeta = tempBeta.multiply(new BigDecimal(i));
			double beta = tempBeta.doubleValue();
			System.out.println("----beta: "+ beta);
			
			long start1 = System.currentTimeMillis();
			PrimitiveMatrix diffusionMatrix = hn2m.createDiffusionMatrixOJA(largestComponent, geneSet, beta);
			long end1 = System.currentTimeMillis();
			System.out.println("\tDiffusion Matrix Time Taken: " + ((end1 - start1) / 1000) + "seconds");

			saveInfluenceGeneCountRange(directory, influenceFile, largestComponent, sourceProteins, diffusionMatrix, geneSet, tempBeta);			
		}
	}

	/**
	 * Saves the results from calculating the betweenness centrality for all nodes in a graph.
	 * @param directory - Directory to save the file in.
	 * @param fileName - File name to save the betweenness centrality results in.
	 * @param graph - Graph used for calculating betweenness centrality.
	 * @throws IOException 
	 */
	private void saveBetweennessCentrality(String directory, String fileName, Graph<String, String> graph) throws IOException{
		File file = new File(directory + "/output/" + fileName);
		file.getParentFile().mkdirs();
		if (!file.exists())
			file.createNewFile();
		String filePath = directory + "/output/" + fileName;
		PrintWriter pw = new PrintWriter(filePath);
		BetweennessCentrality<String, String> bc = new BetweennessCentrality<String, String>(graph);
		bc.setRemoveRankScoresOnFinalize(false);
		bc.evaluate();
		for (String vertex: graph.getVertices())
			pw.println(bc.getVertexRankScore(vertex) + "\t" + vertex);
		pw.close();
	}
	
	/**
	 * Gets 5 quartile source proteins based results in the betweenness centrality file.  
	 * @param directory - Directory where the betweenness centrality score file is located.
	 * @param fileName - Name of the betweenness centrality score file.
	 * @return a set containing the genes from the minimum, 25% quartile, median, 75% quartile, and maximum betweenness centrality. 
	 * @throws IOException
	 */
	private List<String> getSourceProteins(String directory, String fileName) throws IOException{
		File file = new File(directory + "/" + fileName);
		if (!file.exists()){
			System.err.println("betweennessScore.txt file with betweenness scores for all genes in largest network does not exist");
			System.exit(0);
		}

		//read in scores for finding genes in quartile
		Path path = Paths.get(directory, fileName);
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);

		List<Double> bcScores = new ArrayList<Double>(); 
		Map<Double, List<String>> bcMap = new HashMap<Double, List<String>>();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split("\t");
			Double score = Double.parseDouble(row[0]);
			String gene = row[1];
			bcScores.add(score);
			List<String> geneList = bcMap.get(score);
			if (geneList == null){
				geneList = new ArrayList<String>();
				geneList.add(gene);
				bcMap.put(score, geneList);
			}
			else {
				geneList.add(gene);
				bcMap.put(score, geneList);
			}
		}
		br.close();
		Collections.sort(bcScores);
		
		List<Double> quartile = new ArrayList<Double>();
		List<String> quartileGenes = new ArrayList<String>();
		int min = 0;  
		int q25 = (int) Math.ceil(bcScores.size()*0.25);
		int med = (int) Math.ceil(bcScores.size()*0.50);
		int q75 = (int) Math.ceil(bcScores.size()*0.75);
		int max = (int) Math.ceil(bcScores.size()-1);
		quartile.add(bcScores.get(min));
		quartile.add(bcScores.get(q25));
		quartile.add(bcScores.get(med));
		quartile.add(bcScores.get(q75));
		quartile.add(bcScores.get(max));
		for (double q: quartile){
			List<String> geneList = bcMap.get(q);
			Integer size = geneList.size();
			Random randomGen = new Random();
			Integer random = randomGen.nextInt(size);
			quartileGenes.add(geneList.get(random));
		}
		System.out.println("Source Proteins: " + quartileGenes);
		return quartileGenes;
	}

	
	/**
	 * Saves influence gene count.
	 * <p>
	 * Saves the number of proteins with at least a diffusion value of theta influence for the 3 categories below related to the source proteins:
	 * <ul>
	 * <li>direct neighbors (source protein's direct neighbor)
	 * <li>secondary neighbors (source protein's secondary neighbor.  nodes distance 2 away from source protein)
	 * <li>all genes (all genes in network)
	 * </ul>
	 * @param directory - Directory to save files.
	 * @param fileName - File name to save number of genes with at least theta influence for the 3 categories in.
	 * @param graph - Graph used to find neighbors.
	 * @param sourceProteins - List of 5 source proteins obtained from getSourceProteins().
	 * @param diffusionMatrix - Diffusion matrix who's contents are compared with the theta influence value.
	 * @param geneSet- Set of genes that determines matrix ordering.
	 * @param beta - parameter that determines the fraction of own heat each gene retains.
	 * @throws IOException 
	 */
	private void saveInfluenceGeneCountRange(String directory, String fileName, Graph<String,String> graph, List<String> sourceProteins, PrimitiveMatrix diffusionMatrix, Set<String> geneSet, BigDecimal tempBeta) throws IOException{
		Set<String> influenceSet = new TreeSet<String>();
		for (int i=0; i<1001; i++){
			BigDecimal tempPt = new BigDecimal("0.0001");
			tempPt = tempPt.multiply(new BigDecimal(i));
			double influence = tempPt.doubleValue();
			int directNeighborsNum = 0;
			int secondaryNeighborsNum = 0;
			int allGenesNum = 0;
			
			for(String sp: sourceProteins){
				//Obtain the source protein's direct neighbors, secondary neighbors, and all proteins in the network 
				Set<String> directNeighbors = new TreeSet<String>();
				Set<String> secondaryNeighbors = new TreeSet<String>();
				Set<String> allGenes = new TreeSet<String>();
				directNeighbors.addAll(graph.getNeighbors(sp));
				for (String d: directNeighbors)
					secondaryNeighbors.addAll(graph.getNeighbors(d));
				secondaryNeighbors.remove(sp);
				secondaryNeighbors.removeAll(directNeighbors);
				allGenes.addAll(graph.getVertices());
				
				//Get distribution of proteins for the 3 influence point value results
				directNeighborsNum += calculateInfluenceQuantity(diffusionMatrix, geneSet, directNeighbors, sp, influence);
				secondaryNeighborsNum += calculateInfluenceQuantity(diffusionMatrix, geneSet, secondaryNeighbors, sp, influence);
				allGenesNum += calculateInfluenceQuantity(diffusionMatrix, geneSet, allGenes, sp, influence);
				
				//Store the 3 different distribution of proteins
				influenceSet.add(tempBeta +  "\t" + tempPt + "\t" + sp + "\t" + directNeighborsNum + "\t" + secondaryNeighborsNum + "\t" + allGenesNum);
			}
			FileUtils fu = new FileUtils();
			String newDirectory = directory +"/inflectionPoint/"; 
			String newFileName = tempBeta + "_" + fileName;
			fu.saveSetToFile(newDirectory, newFileName, influenceSet);
		}
	}
	
	/**
	 * Calculates the number of proteins greater than or equal to the theta influence point.
	 * @param diffusionMatrix - Diffusion matrix obtained from createDiffusionMatrix(). 
	 * @param geneSet - Set of genes that determines matrix ordering.
	 * @param specificGeneSet - Set of genes who must be found in a matrix.
	 * @param sourceProtein- Source protein used to find matrix value for comparison against the theta influence.
	 * @param thetaInfluence - Parameter used to determine number of genes with at least this value.
	 * @return the number of genes greater than or equal to the theta influence.
	 */
	private int calculateInfluenceQuantity(PrimitiveMatrix diffusionMatrix, Set<String> geneSet, Set<String> specificGeneSet, String sourceProtein, double thetaInfluence){
		int quantity = 0; 
		List<Integer> proteinList = getSpecificGeneIndexInSet(geneSet, specificGeneSet, sourceProtein);
		for (Integer p: proteinList){
			if (proteinList.get(0) != p){
				if (thetaInfluence <= diffusionMatrix.get(proteinList.get(0),p))
					quantity+=1;
			}
		}	
		return quantity; 
	}

	/**
	 * Gets the indices of specific genes within a matrix.
	 * <p>
	 * <b>Note:</b> only used in selectBeta() and matrix is ordered by geneSet.
	 * @param geneSet - Set of genes that determines matrix ordering.
	 * @param specificGeneSet - Set of genes whose indices must be found in a matrix. 
	 * @param sourceProtein - Source protein used to find matrix value for comparison against the theta influence.
	 * @return list of matrix indicies needed for calculateInfleuenceQuantity().
	 */
	private List<Integer> getSpecificGeneIndexInSet(Set<String> geneSet, Set<String> specificGeneSet, String sourceProtein) {
		List<String> geneList = new ArrayList<String>(geneSet);
		List<Integer> indexList = new ArrayList<Integer>();
		indexList.add(geneList.indexOf(sourceProtein));
		for (String s: specificGeneSet){
			indexList.add(geneList.indexOf(s));}
		return indexList;
	}
}