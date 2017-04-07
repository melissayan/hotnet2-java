package edu.ohsu.hotnet2;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

public class HotNet2 {
	
	public static void main(String args[])throws IOException {
		System.out.println("--------------------------- START ---------------------------"); 
		Path currentPath = Paths.get(""); 
		String directory = currentPath.toAbsolutePath().toString();
		String allGenesFile = "allGenes.txt";
		String allPairsFile = "allPairs.txt";
		String largestCompFile = "largestComponentPairs.txt";
//		String refLargestCompFile = "FIsInGene_031516_BigComp_temp.txt";
//		String refLargestCompFile = "FIsInGene_031516_BigComp.txt";
		String betweennessScoreFile = "betweennessScore.txt";
		String influenceFile = "influence.txt";
		String heatScoreFile = "mutation_frequency_expr_filtered.txt";
		String geneMissingHeatFile = "geneMissingHeat.txt";
		String geneIndexFile = "geneIndexReactome.txt";
		String edgeListFile = "edgeListReactome.txt";
		String permuteEdgeListFile = "edgeListPermuted";
		
		BetaSelection bs = new BetaSelection();
		DeltaSelection ds = new DeltaSelection();
		FileUtils fu = new FileUtils();
		GraphUtils gu = new GraphUtils();
		HotNet2Matrix hn2m = new HotNet2Matrix();
		
		if (args.length==0){
			System.err.println("No arguments given");
			System.err.println("Options for [args0]:\n" +
							"\t1\t select beta using “FIsInGene_temp.txt” network\n" +
							"\t2\t select beta using 2015 Reactome FI network\n" +
							"\t\t\tOptional [args1]:\n" +
							"\t\t\tonlyTP53\t only obtain TP53 beta selection results; default will obtain 5 source proteins based on betweeness centrality\n" +
							"\t3\t select beta using “iref_edge_list_temp\n" +
							"\t4\t select beta using network of choice, but file must meet python HotNet2 edgelist requirements\n" +
							"\t\t\tOption for [args1]:\n" +
							"\t\t\t[edgelistFile]\t must provide edgelist file of network\n" +
							"\t5\t select delta using 2015 Reactome FI network\n" +
							"\t\t\tOption for [args1]:\n" +
							"\t\t\t[betaValue]\t must provide a beta parameter, ex. 0.5\n" + 
							"\t6\t select delta using iRefIndex network\n" +
							"\t7\t run HotNet2 without statistical analysis using 2015 Reactome FI network\n" +
							"\t\t\tOption for [args1]:\n" +
							"\t\t\t[betaValue]\t must provide a beta parameter, ex. 0.5\n" + 
							"\t\t\tOption for [args2]:\n" +
							"\t\t\t[deltaValue]\t must provide delta parameter, ex. 0.1\n" +
							"\t8\t generate 2015 Reactome FI network edgelist and gene index for python HotNet2");
			System.exit(0);
		}
		int option = Integer.parseInt(args[0]);
		switch (option) {
		    case 1 : //BetaSelection
		        runBetaSelection(directory, betweennessScoreFile, influenceFile, bs);
		        break;
		    case 2 :
		        runBetaSelection(args, directory, betweennessScoreFile, influenceFile, bs);
		        break;
		    case 3 :
		        runBetaSelection(directory, betweennessScoreFile, influenceFile, bs, fu);
		        break;
		    case 4 :
		    	if (args.length==1){
					System.err.println("Args[1] must indicate the network's edgelist file for beta selection");
					System.err.println("\tNote: edgelist file must meet Python HotNet2's edgelist requirements");
					System.exit(0);
				} 
		    	String networkFile = args[1];
		    	runBetaSelection(directory, networkFile, betweennessScoreFile, influenceFile, bs, fu);
		    	break;
		    case 5 : //DeltaSelection
		        runDeltaSelection(args, directory, ds); 
		        break;
		    case 6 :
		        runDeltaSelection(directory, ds); 
		        break;
		    case 7 : //HotNet2 Algorithm
		        runHotNet2(args, directory, heatScoreFile, hn2m);
		        break;
		    case 8 : //Prepare Files for Python version of HotNet2
		        System.out.println("Preparing Reactome FI network edgelist and gene index for Python version of HotNet2");
	            fu.testFilesForPY();
		        break;
		    default :
		        System.err.println("Don't know the provided option!");
		}
		System.out.println("---------------------------- END ----------------------------");
	}

    private static void runBetaSelection(String directory, String betweennessScoreFile, String influenceFile,
	                                  BetaSelection bs)
	        throws IOException {
    	System.out.println("Beta Selection using FIsInGene_temp.txt network");
	    String fiFile = "FIsInGene_temp.txt";
	    Boolean onlyTP53 = false; 
	    bs.selectBeta(directory, fiFile, betweennessScoreFile, influenceFile, onlyTP53);
	}

	private static void runBetaSelection(String[] args, String directory, String betweennessScoreFile,
	                                  String influenceFile, BetaSelection bs)
	        throws IOException {
	    System.out.println("Beta Selection using Reactome FI network");
	    String fiFile = "FIsInGene_031516_with_annotations.txt";

	    if (args.length==1){
	    	Boolean onlyTP53 = false;
	        bs.selectBeta(directory, fiFile, betweennessScoreFile, influenceFile, onlyTP53);
	    }
	    if (args[1].equals("onlyTP53")){
	        Boolean onlyTP53 = true;
	        bs.selectBeta(directory, fiFile, betweennessScoreFile, influenceFile, onlyTP53);
	    }
	    else {
	        Boolean onlyTP53 = false;
	        bs.selectBeta(directory, fiFile, betweennessScoreFile, influenceFile, onlyTP53);
	    }
	}

	private static void runBetaSelection(String directory, String betweennessScoreFile, String influenceFile,
                                      BetaSelection bs, FileUtils fu)
            throws IOException {
        System.out.println("Beta Selection using iRefIndex network");
        String fiFile = "iref_edge_list";
//        String fiFile = "iref_edge_list_temp";//to use this network, modify BetaSelection.java: comment line 134 and uncomment line 135
        fiFile = fu.modifyEdgeFileToTabDelimited(directory, fiFile, false, " ");
        bs.selectBetaForIrefindex(directory, fiFile, betweennessScoreFile, influenceFile);
    }
	
	private static void runBetaSelection(String directory, String networkFile, String betweennessScoreFile, String influenceFile,
                                      BetaSelection bs, FileUtils fu)
            throws IOException {
		System.out.println("Beta Selection using " + networkFile);
		String fiFile = networkFile;
		fiFile = fu.modifyEdgeFileToTabDelimited(directory, fiFile, false, " ");
		bs.selectBeta(directory, fiFile, betweennessScoreFile, influenceFile, false);
	}

	private static void runDeltaSelection(String[] args, String directory, DeltaSelection ds) throws IOException {
	    System.out.println("Delta Selection using Reactome FI network - 100 permutations");
		if (args.length==1){
			System.err.println("Args[1] must indicate the beta value to be used for delta selection, ex. 0.5");
			System.exit(0);
		}    
	    double beta = Double.parseDouble(args[1]);
	    System.out.println("beta = " + beta);
	    int numPermutations = 100;
	    ds.selectDeltaByCompSize(directory, beta, numPermutations);
	}

	private static void runDeltaSelection(String directory, DeltaSelection ds) throws IOException {
	    System.out.println("Delta Selection using iRefIndex network - 100 permutations");
	    double beta = Double.parseDouble("0.45");
	    System.out.println("beta = " + beta);
	    int numPermutations = 100;
	    ds.selectDeltaForIrefindexByCompSize(directory, beta, numPermutations);
	}

	private static void runHotNet2(String[] args, String directory, String heatScoreFile, HotNet2Matrix hn2m)
	        throws IOException {
	    System.out.println("Reactome FI network - HotNet2 Algorithm without statistical analysis using provided beta and delta parameters - 100 permutations");
	    String fiFile = "FIsInGene_031516_with_annotations.txt";
		if (args.length!=3){
			System.err.println("Args[1] must indicate the beta value, ex. 0.5\n" +
							"Args[2] must indicate the delta value, ex. 0.1");
			System.exit(0);
		}
	    double beta = Double.parseDouble(args[1]);
	    System.out.println("beta = " + beta);
	    double delta = Double.parseDouble(args[2]);
	    System.out.println("delta = " + delta);
	    hn2m.testHotNet2Algorithm(directory, fiFile, heatScoreFile, beta, delta);
	}
}


