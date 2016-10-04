package edu.ohsu.hotnet2;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import org.junit.*;
import org.ojalgo.matrix.PrimitiveMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import edu.uci.ics.jung.algorithms.cluster.WeakComponentClusterer;
import edu.uci.ics.jung.graph.Graph;

public class StatisticalSignificance{
	
	public StatisticalSignificance(){
		
	}
	
	
	private void getStatisticalSignificance(String directory, double beta, double delta, int numPermutations) throws IOException{
		HashMap<Integer, Integer> observed = obtainSubnetworkSizeToCountSumMapReal(directory, beta, delta);
		HashMap<Integer, Integer> expected = obtainSubnetworkSizeToCountSumMapPermuted(directory, beta, delta);
		SortedSet<Integer> sizeSet = new TreeSet<Integer>(observed.keySet());
		sizeSet.addAll(expected.keySet());
		List<Integer> sizeList = new ArrayList<Integer>(sizeSet);
		HashMap<Integer, List<String>> sizeToStatsMap = new HashMap<Integer, List<String>>();
		HashMap<Integer, Double> sizeToPvalueMap = new HashMap<Integer, Double>();
		for (int s: sizeList){
			int observedValue = observed.get(s);
			int expectedValue = expected.get(s);
			double pvalue = observedValue/numPermutations;
			List<String> dataList = new ArrayList<String>();
			dataList.add(Integer.toString(s));
			dataList.add(Double.toString(observedValue));
			dataList.add(Double.toString(expectedValue));
			sizeToStatsMap.put(s, dataList);
			sizeToPvalueMap.put(s, pvalue);
		}
		//save results to file
		String sigDirectory = directory + "/outcome/";
		String sigFile = "significance.txt";
		File file = new File(sigDirectory + sigFile);
		file.getParentFile().mkdirs();
		if (!file.exists())
			file.createNewFile();
		String filePath = sigDirectory + sigFile;
		PrintWriter pw = new PrintWriter(filePath);
		List<Integer> keylist = new ArrayList<Integer>(sizeToStatsMap.keySet()); 
		for(int key: keylist){
			List<String> tempList = sizeToStatsMap.get(key);
			pw.println(key + "\t" + tempList.get(0) + "\t" + tempList.get(1));
		}
		pw.close();
	}
	
	private HashMap<Integer, Integer> obtainSubnetworkSizeToCountSumMapReal(String directory, double beta, double delta) throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		HashMap<Integer, Integer> sumSizeCountMap = new HashMap<Integer, Integer>();
		//Create set of genes in network and HashMap of gene index to allow gene names in graph instead of numbers
		String geneIndexDirectory = directory;
		String geneIndexFile = "geneIndexReactome.txt";
		Set<String> genes = fu.getAllGenesPY(geneIndexDirectory, geneIndexFile, "\t");
		HashMap<String, String> geneIndexMap = fu.getGeneIndexNamePY(geneIndexDirectory, geneIndexFile, "\t");
		//Read in all permutation files within a directory to get subnetwork size sum count
		HashMap<Integer, Integer> sizeToCountMap = new HashMap<Integer, Integer>();
		String edgeListDirectory = directory;
		File folder = new File(edgeListDirectory);
		File[] listOfFiles = folder.listFiles();
		for (File file : folder.listFiles()){
			String edgeListFile = file.getName();
			String heatScoreDirectory = directory;
			String heatScoreFile = "heatScore.txt";
			//Create graph from file and check there is only 1 component
			Set<String> pairs = fu.getAllInteractionPairsPY(edgeListDirectory, edgeListFile, geneIndexMap, "\t");
			Graph<String, String> largestComponent = gu.createGraph(genes, pairs);
			WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
			Set<Set<String>> components = wcc.transform(largestComponent);
			if (components.size() != 1){
				throw new IllegalArgumentException("Provided permuted graph " + edgeListFile + " is not connected, it has " + components + " components.");
			}
			//Get sum of subnetwork counts from individual networks
			HashMap<Integer, Integer> indvidualSizeCount = obtainIndividualSubnetworkSizeCount(heatScoreDirectory, heatScoreFile, largestComponent, beta, delta);
			List<Integer> sizeList = new ArrayList<Integer>(sizeToCountMap.keySet());
			for (int s: sizeList){
				int count = sizeToCountMap.get(s);
				if (sizeToCountMap.get(s) == null)
					sizeToCountMap.put(s, 1);
				else
					sizeToCountMap.put(s, count+indvidualSizeCount.get(s));
			}
		}
		//Save subnetwork size sum count
		String sigDirectory = directory + "/outcome/";
		String sigFile = "significance.txt";
		fu.saveHashMapToFile(sigDirectory, sigFile, sizeToCountMap);
		return sizeToCountMap;
	}
	
	//For all random networks, get a sum of greater or equal counts of subnetwork sizes
	private HashMap<Integer, Integer> obtainSubnetworkSizeToCountSumMapPermuted(String directory, double beta, double delta) throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		HashMap<Integer, Integer> sumSizeCountMap = new HashMap<Integer, Integer>();
		//Create set of genes in network and HashMap of gene index to allow gene names in graph instead of numbers
		String geneIndexDirectory = directory + "/PythonPermutation/";
		String geneIndexFile = "geneIndexReactome.txt";
		Set<String> genes = fu.getAllGenesPY(geneIndexDirectory, geneIndexFile, "\t");
		HashMap<String, String> geneIndexMap = fu.getGeneIndexNamePY(geneIndexDirectory, geneIndexFile, "\t");
		//Read in all permutation files within a directory to get subnetwork size sum count
		HashMap<Integer, Integer> sizeToCountMap = new HashMap<Integer, Integer>();
		String edgeListDirectory = directory + "/PythonPermutation/permutations/";
		File folder = new File(edgeListDirectory);
		File[] listOfFiles = folder.listFiles();
		for (File file : folder.listFiles()){
			String edgeListFile = file.getName();
			String heatScoreDirectory = directory+"";
			String heatScoreFile = "heatScore.txt";
			//Create graph from file and check there is only 1 component
			Set<String> pairs = fu.getAllInteractionPairsPY(edgeListDirectory, edgeListFile, geneIndexMap, "\t");
			Graph<String, String> largestComponent = gu.createGraph(genes, pairs);
			WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
			Set<Set<String>> components = wcc.transform(largestComponent);
			if (components.size() != 1){
				throw new IllegalArgumentException("Provided permuted graph " + edgeListFile + " is not connected, it has " + components + " components.");
			}
			//Get sum of subnetwork counts from individual networks
			HashMap<Integer, Integer> indvidualSizeCount = obtainIndividualSubnetworkSizeCount(heatScoreDirectory, heatScoreFile, largestComponent, beta, delta);
			List<Integer> sizeList = new ArrayList<Integer>(sizeToCountMap.keySet());
			for (int s: sizeList){
				int count = sizeToCountMap.get(s);
				if (sizeToCountMap.get(s) == null)
					sizeToCountMap.put(s, 1);
				else
					sizeToCountMap.put(s, count+indvidualSizeCount.get(s));
			}
		}
		//Save subnetwork size sum count
		String sigDirectory = directory + "/outcome/";
		String sigFile = "significance.txt";
		fu.saveHashMapToFile(sigDirectory, sigFile, sizeToCountMap);
		return sizeToCountMap;
	}
	
	
	//For each network, get the counts of subnetwork sizes
	/**
	 * Obtains the number of each subnetwork size found in the provided graph. 
	 * @param heatScoreDirectory - Directory of heat score file location. 
	 * @param heatScoreFile - File name with no header containing genes and heat scores separated a space. 
	 * @param graph - Graph used to extract subnetworks from. 
	 * @param beta - Fraction of own heat each gene retains. 
	 * @param delta - Minimum edge weight threshold; values below edgeWeight are set to 0.
	 * @return a HashMap with the subnetwork size as key and the number of it's occurrences as the value. 
	 * @throws IOException
	 */
	private HashMap<Integer, Integer> obtainIndividualSubnetworkSizeCount(String heatScoreDirectory, String heatScoreFile, Graph<String, String> graph, double beta, double delta) throws IOException{
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		HotNet2Matrix hn2m = new HotNet2Matrix();
		Set<String> geneSet = gu.getGeneGraphSet(graph);
		PrimitiveMatrix tempF = hn2m.createDiffusionMatrixOJA(graph, geneSet, beta);
		RealMatrix F = hn2m.convertOJAToACM(tempF);	
		RealMatrix E = hn2m.createExchangedHeatMatrix(heatScoreDirectory, heatScoreFile, F, geneSet);
		RealMatrix H = hn2m.identifyHotSubnetworks(E, delta);
		Graph<String, String> graphFromMatrix = gu.covertMatrixToDirectedGraph(H, geneSet);
		WeakComponentClusterer<String, String> wcc = new WeakComponentClusterer<String, String>();		
		Set<Set<String>> components = wcc.transform(graphFromMatrix);
		HashMap<Integer, Integer> sizeToCountMap = new HashMap<Integer, Integer>();
		for (Set<String> c: components){
			int size = c.size();
			int count = sizeToCountMap.get(size);
			if (sizeToCountMap.get(size) == null)
				sizeToCountMap.put(size, 1);
			else
				sizeToCountMap.put(size, count+1);
		}
		return sizeToCountMap;
	}

}