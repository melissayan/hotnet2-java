package edu.ohsu.hotnet2;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.DoubleStream;

import org.junit.*;

import edu.uci.ics.jung.graph.*;//Graph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.algorithms.cluster.*;//WeakComponentClusterer;
import edu.uci.ics.jung.algorithms.filters.*;//FilterUtils;
import edu.uci.ics.jung.algorithms.importance.*;//BetweennessCentrality;
import edu.uci.ics.jung.algorithms.matrix.*;
import edu.uci.ics.jung.graph.*;

import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SparseFieldMatrix;
import org.apache.commons.math3.stat.*;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.UnivariateStatistic;
import org.ojalgo.*;
import org.ojalgo.OjAlgoUtils;
import org.ojalgo.access.Access2D.Builder;
import org.ojalgo.matrix.*;
import org.ojalgo.matrix.jama.*;
import org.ojalgo.matrix.store.PrimitiveDenseStore;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;

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
		
		if (args.length==0){
			System.err.println("Arg[0] must indicate one of the ReactomeFI files below using a number:");
			System.err.println("\t 1: FIsInGene_temp.txt\n\t 2: FIsInGene_031516_with_annotations.txt");
			System.exit(0);
		}
		if (Integer.parseInt(args[0])==1){
			String fiFile = "FIsInGene_temp.txt";
			Path fiFilePath = Paths.get(directory, fiFile);
			bs.selectBeta(directory, fiFile, betweennessScoreFile, influenceFile);
		}
		if (Integer.parseInt(args[0])==2){
			String fiFile = "FIsInGene_031516_with_annotations.txt";
			Path fiFilePath = Paths.get(directory, fiFile);
//			bs.selectBeta(directory, fiFile, betweennessScoreFile, influenceFile);
			Graph<String, String> allGenesGraph = gu.createReactomeFIGraph(directory, fiFile);
			Graph<String, String> largestComponent = gu.createLargestComponentGraph(allGenesGraph);
			System.out.println("size: " + largestComponent.getVertices().size());
		}
		if (Integer.parseInt(args[0])==3){
			String fiFile = "iref_edge_list";
//			String fiFile = "iref_edge_list_temp";
			fiFile = fu.modifyEdgeFileToTabDelimited(directory, fiFile, false, " ");
			System.out.println(fiFile);
			Path fiFilePath = Paths.get(directory, fiFile);
			bs.selectBetaForIrefindex(directory, fiFile, betweennessScoreFile, influenceFile);
		}
		//DeltaSelection
		if (Integer.parseInt(args[0])==4){
			System.out.println("Reactome FI network");
			double beta = Double.parseDouble(args[1]);
			int numPermutations = 100;
			ds.selectDeltaByCompSize(directory, beta, numPermutations);			
		}
		if (Integer.parseInt(args[0])==5){
			System.out.println("iRefIndex network");
			double beta = Double.parseDouble("0.45");
			int numPermutations = 100;
			ds.selectDeltaForIrefindexByCompSize(directory, beta, numPermutations);			
		}
		
		System.out.println("---------------------------- END ----------------------------");
	}
}

