package hotnet2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.junit.Test;

import edu.uci.ics.jung.graph.Graph;

public class FileUtils{

	public FileUtils(){
		
	}
	
	@Test
	public void testFilesForPY() throws IOException {
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
		String largestCompFile = "largestComponentPairs.txt";
		String geneIndexFile = "geneIndexReactome.txt";
		String edgeListFile = "edgeListReactome.txt";
				
		convertFileForPythonHotNet2(directory, largestCompFile, geneIndexFile, edgeListFile);
	}
	
	@Test
	public void testGetAllGenesReactome() throws IOException{
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
//		String fiFile = "FIsInGene_031516_with_annotations.txt";
		String fiFile = "FIsInGene_031516_BigComp.txt";
		
		Path fiFilePath = Paths.get(directory, fiFile);
		Set<String> allGenes = getAllGenesReactome(fiFilePath);
		System.out.println(allGenes.size());
	}
	
	@Test
	public void testGetPYGenes() throws IOException {
		String fileName = "geneIndexReactome.txt";
		FileReader fileReader = new FileReader(fileName);
		BufferedReader br = new BufferedReader(fileReader);
		String line = null;
		Set<String> genes = new HashSet<String>();
		while ((line = br.readLine()) != null) {
			String[] tokens = line.split("\t");
			if (tokens[1].contains(" "))
				continue;
			genes.add(tokens[1]);
		}
		br.close();
		System.out.println("Total genes: " + genes.size());
	}
	
	@Test
	public void testCompare2GeneIndexFiles() throws IOException{
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
		String file1 = "geneIndexReactome-old20160824.txt"; //JavaResults
		String file2 = "reactome_index_genes"; //PythonResults
		
		//Store genes from file1 and file2 into separate sets
		Path path1 = Paths.get(directory, file1);
		Set<String> geneFile1 = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path1, charset);
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split("\t");
			String gene = row[1];
			geneFile1.add(gene);
		}
		br.close();
		Path path2 = Paths.get(directory, file2);
		Set<String> geneFile2 = new HashSet<String>();
		br = Files.newBufferedReader(path2, charset);
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(" ");
			String gene = row[1];
			geneFile2.add(gene);
		}
		br.close();
		
		//Find differences
		Set<String> temp1 = new HashSet<String>(geneFile1);
		Set<String> temp2 = new HashSet<String>(geneFile2);
		geneFile1.removeAll(temp2);
		System.out.println("Number of Genes in " + file1 + " and not in " + file2 + ": " + geneFile1.size());
		System.out.println("\t" + geneFile1);
		
		geneFile2.removeAll(temp1);
		System.out.println("Number of Genes in " + file2 + " and not in " + file1 + ": " + geneFile2.size());
		System.out.println("\t" + geneFile2);
	}
	
	/**
	 * Gets file contents and stores into a set.
	 * @param directory - Directory containing file.	
	 * @param fileName - Name of file to read in.
	 * @return a ordered set of file line contents.
	 * @throws IOException
	 */
	public Set<String> getFileContents (String directory, String fileName) throws IOException{
		Set<String> contents = new TreeSet<String>();
		Path path = Paths.get(directory, fileName);
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		for (String line = null; (line = br.readLine()) != null;){
			contents.add(line);
		}
		br.close();
		return contents;
	}
	
	/**
	 * Modifies a file containing gene pairs by saving it as a new tab delimited file.
	 * <p>
	 * <b>Note:</b> used to test out interactions from the Python implementation.
	 * @param directory	 - Directory containing file to be modified.
	 * @param fileToModify - File that will be modified into a tab delimited file.
	 * @param header - Boolean to determine if the file has a header (true) or no header (false).
	 * @param delimiter - Delimiter in the file that will be modified.
	 * @return the name of the newly modified file.
	 * @throws IOException
	 */
	private String modifyEdgeFileToTabDelimited(String directory, String fileToModify, Boolean header, String delimiter) throws IOException{
		//Set up for file that will get modified
		Path path = Paths.get(directory, fileToModify);
		Set<String> genes = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		if (header == true)
			br.readLine();
		//Set up for modified file
		String modifiedFileName = "modified_" + fileToModify; 
		File file = new File(modifiedFileName);		
		if (!file.exists())
			file.createNewFile();
		String filePath = directory +"/" + modifiedFileName;
		PrintWriter pw = new PrintWriter(filePath);
		pw.println("gene1 \tgene2");
		//Read file that 
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(delimiter);
			String gene1 = row[0];
			String gene2 = row[1];
			pw.println(gene1 + "\t" + gene2);
		}
		br.close();
		pw.close();
		return modifiedFileName;
	}
	
	/**
	 * Saves a set into a file.
	 * @param directory - Directory to save the file in.	
	 * @param fileName - Name of file to save set in.
	 * @param setSave - Set to be saved. 
	 * @throws IOException 
	 */
	public void saveSetToFile(String directory, String fileName, Set<String> setSave) throws IOException{
		File file = new File(directory + fileName);
		file.getParentFile().mkdirs();
		if (!file.exists())
			file.createNewFile();
		String filePath = directory + fileName;
		PrintWriter pw = new PrintWriter(filePath);
		for(String s: setSave){
			pw.println(s);
		}
		pw.close();
	}
	
	/**
	 * Compares 2 different files for similarities and differences.
	 * <p>
	 * <b>Note:</b> files must be in same directory. 
	 * @param directory - Directory containing files
	 * @param file1 - Name of the first file to read in
	 * @param file2 - Name of the second file to read in
	 * @throws IOException
	 */
	public void compareFiles (String directory, String file1, String file2) throws IOException{
		Set<String> content1 = getFileContents(directory, file1);
		Set<String> content2 = getFileContents(directory, file2);
		if (content1.equals(content2) == true)
			System.out.println(file1 + " and " + file2 + "have the same contents");
		else{
			System.out.println("The following are in " + file1 + " but not in " + file2 +":");
			Set<String> temp = content1; 
			temp.removeAll(content2);
			Integer i = 1;
			for (String t: temp){
				System.out.println(i + ".\t" + t);
				i++; 
			}	
		}
	}

	
	/**
	 * Gets all genes from the Reactome FI network text file.
	 * @param path - Path to file location.
	 * @return a set containing all genes
	 * @throws IOException
	 */
	public Set<String> getAllGenesReactomeWrapper(Path path) throws IOException{
		return getAllGenesReactome(path);
	}
	
	/**
	 * Gets all genes from the Reactome FI network text file.
	 * @param path - Path to file location.
	 * @return a set containing all genes
	 * @throws IOException
	 */
	private Set<String> getAllGenesReactome(Path path) throws IOException{
		Set<String> genes = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		br.readLine();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split("\t");
			String gene1 = row[0];
			genes.add(gene1);
			String gene2 = row[1];
			genes.add(gene2);
		}
		br.close();
		return genes;
	}
	
	/**
	 * Gets all interaction pairs from the Reactome FI network text file.
	 * @param path - Path to the file location.
	 * @return a set containing all interaction pairs.
	 * @throws IOException
	 */
	public Set<String> getAllInteractionPairsReactomeWrapper(Path path) throws IOException{
		return getAllInteractionPairsReactome(path);
	}
	
	/**
	 * Gets all interaction pairs from the Reactome FI network text file.
	 * @param path - Path to the file location.
	 * @return a set containing all interaction pairs.
	 * @throws IOException
	 */
	private Set<String> getAllInteractionPairsReactome(Path path) throws IOException{
		Set<String> pairs = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		br.readLine();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split("\t");
			String gene1 = row[0];
			String gene2 = row[1];
			pairs.add(gene1 + "\t"+ gene2);
		}
		br.close();
		return pairs;
	}

	/**
	 * Gets all genes from any network text file with specified delimiter and header.
	 * @param path - Path to file location.
	 * @param delimiter - Delimiter that separates columns in file.
	 * @return a set containing all genes.
	 * @throws IOException
	 */
	public Set<String> getAllGenesWrapper(Path path, String delimiter) throws IOException{
		return getAllGenes(path, delimiter);
	}
	
	/**
	 * Gets all genes from any network text file with specified delimiter and header.
	 * @param path - Path to file location.
	 * @param delimiter - Delimiter that separates columns in file.
	 * @return a set containing all genes.
	 * @throws IOException
	 */
	private Set<String> getAllGenes(Path path, String delimiter) throws IOException{
		Set<String> genes = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		br.readLine();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(delimiter);
			String gene1 = row[0];
			genes.add(gene1);
			String gene2 = row[1];
			genes.add(gene2);
		}
		br.close();
		return genes;
	}

	/**
	 * Gets all interaction pairs from any network text file with specified delimiter and header. 
	 * @param path - Path to the file location.
	 * @param delimiter - Delimiter that separates columns in file.
	 * @return a set containing all interaction pairs.
	 * @throws IOException
	 */
	public Set<String> getAllInteractionPairsWrapper(Path path, String delimiter) throws IOException{
		return getAllInteractionPairs(path, delimiter);
	}
	
	/**
	 * Gets all interaction pairs from any network text file with specified delimiter and header. 
	 * @param path - Path to the file location.
	 * @param delimiter - Delimiter that separates columns in file.
	 * @return a set containing all interaction pairs.
	 * @throws IOException
	 */
	private Set<String> getAllInteractionPairs(Path path, String delimiter) throws IOException{
		Set<String> pairs = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		br.readLine();
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split(delimiter);
			String gene1 = row[0];
			String gene2 = row[1];
			pairs.add(gene1 + "\t"+ gene2);
		}
		br.close();
		return pairs;
	}

	/**
	 * Gets all genes from network text file without a header.
	 * @param path - Path to file location.
	 * @return a set containing all genes.
	 * @throws IOException
	 */
	public Set<String> getAllGenesWrapper(Path path) throws IOException{
		return getAllGenes(path);
	}
	
	/**
	 * Gets all genes from network text file without a header.
	 * @param path - Path to file location.
	 * @return a set containing all genes.
	 * @throws IOException
	 */
	private Set<String> getAllGenes(Path path) throws IOException{
		Set<String> genes = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split("\t");
			String gene1 = row[0];
			genes.add(gene1);
			String gene2 = row[1];
			genes.add(gene2);
		}
		br.close();
		return genes;
	}

	/**
	 * Gets all interaction pairs from network text file without a header.
	 * @param path - Path to the file location.
	 * @return a set containing all interaction pairs.
	 * @throws IOException
	 */
	public Set<String> getAllInteractionPairsWrapper(Path path) throws IOException{
		return getAllInteractionPairs(path);
	}
	
	/**
	 * Gets all interaction pairs from network text file without a header.
	 * @param path - Path to the file location.
	 * @return a set containing all interaction pairs.
	 * @throws IOException
	 */
	private Set<String> getAllInteractionPairs(Path path) throws IOException{
		Set<String> pairs = new HashSet<String>();
		Charset charset = Charset.forName("UTF-8");
		BufferedReader br = Files.newBufferedReader(path, charset);
		for (String line = null; (line = br.readLine()) != null;){
			String[] row = line.split("\t");
			String gene1 = row[0];
			String gene2 = row[1];
			pairs.add(gene1 + "\t"+ gene2);
		}
		br.close();
		return pairs;
	}

	/**
	 * Converts the largest component file into files suitable for running the Python version of HotNet2.
	 * <p>
	 * <b>Note:</b> HotNet2 Python version requires a edge list file and a gene index file.  The HotNet2 CreatePPRMat.py run() method skips whitespace so only part of gene name is kept, to ensure Python code doesn't need changing gene names with spaces are replaced with "_".
	 * @param directory	- Directory where the largest component of the graph is stored as a file of gene pairs.
	 * @param largestCompFile - File name containing interaction pairs from the largest component.
	 * @param geneIndexFile - File name where the index and gene will be saved.
	 * @param edgeListFile - File name where the gene pairs will be saved.
	 * @throws IOException
	 */
	private void convertFileForPythonHotNet2(String directory, String largestCompFile, String geneIndexFile, String edgeListFile) throws IOException{
		Path largeCompFilePath = Paths.get(directory+"/output/", largestCompFile);

		//Get all genes and all interaction pairs from largest Component file
		Set<String> allGenes = getAllGenes(largeCompFilePath);
		Set<String> allPairs = getAllInteractionPairs(largeCompFilePath);	
		
		//Create graph and get the largest component (this checks to ensure file only has 1 component)
		GraphUtils gu = new GraphUtils();
		Graph<String, String> allGenesGraph = gu.createGraphWrapper(allGenes, allPairs);
		Graph<String, String> largestComponent = gu.createLargestComponentGraphWrapper(allGenesGraph);
		Set<String> geneSet = gu.getGeneGraphSetWrapper(largestComponent);
		Set<String> genePairs = gu.getPairsWrapper(largestComponent);
		
		//Create gene index file and a hash map with the gene and index for use in creating the edge list file
		File geneFile = new File(directory + "/output/" + geneIndexFile);
		geneFile.getParentFile().mkdirs();
		if (!geneFile.exists())
			geneFile.createNewFile();
		String geneFilePath = directory + "/output/" + geneIndexFile;
		PrintWriter gPw = new PrintWriter(geneFilePath);
		int i=1;
		Map<String, Integer> geneMap = new HashMap<String, Integer>();
		for(String gs: geneSet){
			String gsModified = gs.replace(" ", "_");	//Python HotNet2 CreatePPRMat.py run() method skips whitespaces, so only 1st part of gene name is kept 
			gPw.println(i + "\t" + gsModified);
			geneMap.put(gs, i);
			i++;
		}
		gPw.close();
		
		//Create edge list file
		File edgeFile = new File(directory + "/output/" + edgeListFile);
		edgeFile.getParentFile().mkdirs();
		if (!edgeFile.exists())
			edgeFile.createNewFile();
		String edgeFilePath = directory + "/output/" + edgeListFile;
		PrintWriter ePw = new PrintWriter(edgeFilePath);
		for(String gp: genePairs){
			String[] pair = gp.split("\t");
			String gene1 = pair[0];
			String gene2 = pair[1];
			ePw.println(geneMap.get(gene1) + "\t" + geneMap.get(gene2));
		}
		ePw.close();
	}
	
}