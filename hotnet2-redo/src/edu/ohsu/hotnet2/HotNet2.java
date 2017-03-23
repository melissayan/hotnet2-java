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
		String matrixOrderFile = "matrixOrder.txt";
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
			System.err.println("Arg[0] must indicate one of the ReactomeFI files below using a number:");
			System.err.println("\t 1: FIsInGene_temp.txt\n\t 2: FIsInGene_031516_with_annotations.txt");
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
		    case 4 : //DeltaSelection
		        runDeltaSelection(args, directory, ds); 
		        break;
		    case 5 :
		        runDeltaSelection(directory, ds); 
		        break;
		    case 6 : //HotNet2 Algorithm
		        runHotNet2(args, directory, heatScoreFile, hn2m);
		        break;
		    case 7 : //Prepare Files for Python version of HotNet2
		        System.out.println("Preparing Reactome FI network edgelist and gene index for Python version of HotNet2");
	            fu.testFilesForPY();
		        break;
		    default :
		        System.err.println("Don't know the provided option!");
		}
		System.out.println("---------------------------- END ----------------------------");
	}

    private static void runHotNet2(String[] args, String directory, String heatScoreFile, HotNet2Matrix hn2m)
            throws IOException {
        System.out.println("Reactome FI network - HotNet2 Algorithm using provided beta and delta parameters");
        String fiFile = "FIsInGene_031516_with_annotations.txt";
        double beta = Double.parseDouble(args[1]);
        double delta = Double.parseDouble(args[2]);
        hn2m.testHotNet2Algorithm(directory, fiFile, heatScoreFile, beta, delta);
    }

    private static void runDeltaSelection(String directory, DeltaSelection ds) throws IOException {
        System.out.println("iRefIndex network");
        double beta = Double.parseDouble("0.45");
        int numPermutations = 100;
        ds.selectDeltaForIrefindexByCompSize(directory, beta, numPermutations);
    }

    private static void runDeltaSelection(String[] args, String directory, DeltaSelection ds) throws IOException {
        System.out.println("Reactome FI network");
        double beta = Double.parseDouble(args[1]);
        System.out.println("beta = " + beta);
        int numPermutations = 100;
        ds.selectDeltaByCompSize(directory, beta, numPermutations);
    }

    private static void runBetaSelection(String directory, String betweennessScoreFile, String influenceFile,
                                      BetaSelection bs, FileUtils fu)
            throws IOException {
        System.out.println("Beta Selection using iRefIndex network");
        String fiFile = "iref_edge_list";
//	          String fiFile = "iref_edge_list_temp";
        fiFile = fu.modifyEdgeFileToTabDelimited(directory, fiFile, false, " ");
        System.out.println(fiFile);
        bs.selectBetaForIrefindex(directory, fiFile, betweennessScoreFile, influenceFile);
    }

    private static void runBetaSelection(String[] args, String directory, String betweennessScoreFile,
                                      String influenceFile, BetaSelection bs)
            throws IOException {
        System.out.println("Beta Selection using Reactome FI network");
        String fiFile = "FIsInGene_031516_with_annotations.txt";
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
                                      BetaSelection bs)
            throws IOException {
        String fiFile = "FIsInGene_temp.txt";
        Boolean onlyTP53 = false; 
        bs.selectBeta(directory, fiFile, betweennessScoreFile, influenceFile, onlyTP53);
    }
}


