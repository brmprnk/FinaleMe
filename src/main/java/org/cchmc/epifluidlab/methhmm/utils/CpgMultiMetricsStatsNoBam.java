/**
 * CpgMultiMetricsStats.java
 * Feb 27, 2016
 * 5:28:37 PM
 * yaping    lyping1986@gmail.com
 */
package org.cchmc.epifluidlab.methhmm.utils;


import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.IntervalTree.Node;
import htsjdk.tribble.CloseableTribbleIterator;
import htsjdk.tribble.TabixFeatureReader;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.Pair;
import org.apache.log4j.Logger;
import org.biojava.nbio.genome.parsers.twobit.TwoBitParser;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import edu.unc.genomics.io.BigWigFileReader;

/**
 *
 */
public class CpgMultiMetricsStatsNoBam {

	@Option(name="-kmerLen",usage="the K-mer length to check. Default: 4")
	public int kmerLen = 4;

	@Option(name="-kmerString",usage="the fiel contain selected K-mer to check. Otherwise, use -kmerLen to automately generate all the k-mer. Default: null")
	public String kmerString = null;

	@Option(name="-kmerExt",usage="the +/- region in reference genome to check the k-mer frequency. default is +/- 100bp around CpGs. Default: 100")
	public int kmerExt = 100;
	
	@Option(name="-excludeRegions",usage="bed files indicated excluded regions. -excludeRegions trackFileName. Default: null")
	public ArrayList<String> excludeRegions = null;

	@Option(name="-overlapRegions",usage="bed files to check if regions are overlapped. -overlapRegions trasckName:trackFileName. Default: null")
	public ArrayList<String> overlapRegions = null;

	@Option(name="-distantRegions",usage="bed files to check the distance to these regions. like the distance to TSS. -distantRegions trasckName:trackFileName. Default: null")
	public ArrayList<String> distantRegions = null;

	@Option(name="-valueWigs",usage="bigwig files to check the value in these regions. like the recombination rate of some region.extRegion indicate the +/- regionBp from CpG, -valueWigs trasckName:extRegion:trackFileName. Default: null")
	public ArrayList<String> valueWigs = null;

	@Option(name="-valueBeds",usage="tabixed bed.gz files to check the value in these regions. like the recombination rate of some region. -valueBeds trasckName:extRegion:trackFileName. Default: null")
	public ArrayList<String> valueBeds = null;

	@Option(name="-skipCpgDist",usage="Skip the calcualtion of CpG distance. Default: false")
	public boolean skipCpgDist = false;

	@Option(name="-h",usage="show option information")
	public boolean help = false;

	
	@Argument
	private List<String> arguments = new ArrayList<String>();

	final private static String USAGE = "CpgMultiMetricsStats [opts] hg19.2bit cpg_list.bed all_cpg.bed cpg_detail.txt.gz";
	
	private static Logger log = Logger.getLogger(CpgMultiMetricsStatsNoBam.class);

	private static long startTime = -1;
	private static long points = 0;


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		CpgMultiMetricsStatsNoBam cmms = new CpgMultiMetricsStatsNoBam();
		//BasicConfigurator.configure();
		cmms.doMain(args);
	}
	
	@SuppressWarnings("resource")
	public void doMain(String[] args)
			throws Exception {

					CmdLineParser parser = new CmdLineParser(this);
					//parser.setUsageWidth(80);
					try
					{
						if(help || args.length < 4) throw new CmdLineException(parser, USAGE, new Throwable());
						parser.parseArgument(args);
						
					
					}
					catch (CmdLineException e)
					{
						System.err.println(e.getMessage());
						// print the list of available options
						parser.printUsage(System.err);
						System.err.println();
						return;
					}

					//read input bed file, for each row,
					//String intervalFile = arguments.get(0);
					String refFile = arguments.get(0);
					String cpgListFile = arguments.get(1);
					String allCpgFile = arguments.get(2);
					String detailFile = arguments.get(3);

					initiate();			
					
					//initiate different kinds of reader
					//reference genome
					TwoBitParser refParser = new TwoBitParser(new File(refFile));
					
					//String[] names = p.getSequenceNames();
					//for(int i=0;i<names.length;i++) {
					//  p.setCurrentSequence(names[i]);
					//  p.printFastaSequence();
					//  p.close();
					//}
					//loading exlusion interval file

					//load interval files
					log.info("Processing interval file ... ");
					
					
					HashMap<String,IntervalTree<Integer>> ignoreLocCollections = null;
					if(excludeRegions!=null && !excludeRegions.isEmpty()){
						log.info("Excluding intervals ... ");
						ignoreLocCollections = new HashMap<String,IntervalTree<Integer>>();
						
						for(String excludeRegion : excludeRegions){
							BufferedReader br = new BufferedReader(new FileReader(excludeRegion));
							String line;
							
							while( (line = br.readLine()) != null){
								if(line.startsWith("#"))
									continue;
								String[] splitLines = line.split("\t");
								if(splitLines.length<3){
									continue;
								}
								String chr = splitLines[0];
								int start = Integer.parseInt(splitLines[1]);
								int end = Integer.parseInt(splitLines[2]);
								IntervalTree<Integer> tree;
								if(ignoreLocCollections.containsKey(chr)){
									tree = ignoreLocCollections.get(chr);
								}else{
									tree = new IntervalTree<Integer>();
								}
								tree.put(start, end, 1);
								ignoreLocCollections.put(chr, tree);
							}
							br.close();
						
						}
					}
					
					
						HashMap<String,IntervalTree<String>> allCpgLocCollections = new HashMap<String,IntervalTree<String>>();
						if(!skipCpgDist){
							log.info("Loading all CpG intervals ... ");
							GZIPInputStream gzipInputStream = null;
							BufferedReader br1;
							if(allCpgFile.endsWith(".gz")){
								gzipInputStream = new GZIPInputStream(new FileInputStream(allCpgFile));
								br1 = new BufferedReader(new InputStreamReader(gzipInputStream));
								
							}else{
								br1 = new BufferedReader(new FileReader(allCpgFile));
							}
								
								String line1;
								
								while( (line1 = br1.readLine()) != null){
									if(line1.startsWith("#"))
										continue;
									String[] splitLines = line1.split("\t");
									if(splitLines.length<3){
										continue;
									}
									String chr = splitLines[0];
									int start = Integer.parseInt(splitLines[1]);
									int end = Integer.parseInt(splitLines[2]);
									IntervalTree<String> tree;
									
									if(allCpgLocCollections.containsKey(chr)){
										tree = allCpgLocCollections.get(chr);
									}else{
										tree = new IntervalTree<String>();
									}
									String strand = ".";
									if(splitLines.length >= 6){
										if(splitLines[5].equalsIgnoreCase("-")){
											strand = "-";
										}else if(splitLines[5].equalsIgnoreCase("+")){
											strand = "+";
										}
										//strand = splitLines[5].equalsIgnoreCase("-") ? "-" : "+";
									}
									tree.put(start, end, strand);
									allCpgLocCollections.put(chr, tree);
								}
								if(allCpgFile.endsWith(".gz")){
									gzipInputStream.close();
								}
								br1.close();
						}
						
						
						
						
					HashMap<String, HashMap<String,IntervalTree<Integer>>> overlapLocStringCollections = null;
					LinkedHashSet<String> overlapLocString = new LinkedHashSet<String>();
					if(overlapRegions!=null && !overlapRegions.isEmpty()){
						log.info("Overalpped intervals ... ");
						overlapLocStringCollections = new HashMap<String, HashMap<String,IntervalTree<Integer>>>();
						for(String overlapRegionString : overlapRegions){
							
							String[] splitStrings = overlapRegionString.split(":");
							if(splitStrings.length < 2){
								throw new IllegalArgumentException("need to provide trackname:trackFile for overlapRegions");
							}
							HashMap<String,IntervalTree<Integer>> overlapLocCollections =  new HashMap<String,IntervalTree<Integer>>();
							String overlapRegionName = splitStrings[0];
							String overlapRegion = splitStrings[1];
							overlapLocString.add(overlapRegionName);
							BufferedReader br = new BufferedReader(new FileReader(overlapRegion));
							String line;
							while( (line = br.readLine()) != null){
								if(line.startsWith("#"))
									continue;
								String[] splitLines = line.split("\t");
								if(splitLines.length<3){
									continue;
								}
								String chr = splitLines[0];
								int start = Integer.parseInt(splitLines[1]);
								int end = Integer.parseInt(splitLines[2]);
								IntervalTree<Integer> tree;
								if(overlapLocCollections.containsKey(chr)){
									tree = overlapLocCollections.get(chr);
								}else{
									tree = new IntervalTree<Integer>();
								}
								tree.put(start, end, 1);
								overlapLocCollections.put(chr, tree);
								
							}
							br.close();
							overlapLocStringCollections.put(overlapRegionName, overlapLocCollections);
						
						}
					}
					
					HashMap<String, HashMap<String,IntervalTree<String>>> distantLocStringCollections = null;
					LinkedHashSet<String> distantLocString = new LinkedHashSet<String>();
					if(distantRegions!=null && !distantRegions.isEmpty()){
						log.info(" Intervals used to calculate distances... ");
						distantLocStringCollections = new HashMap<String, HashMap<String,IntervalTree<String>>>();
						for(String distantRegionString : distantRegions){
							
							String[] splitStrings = distantRegionString.split(":");
							if(splitStrings.length < 2){
								throw new IllegalArgumentException("need to provide trackname:trackFile for distantRegions");
							}
							HashMap<String,IntervalTree<String>> distantLocCollections =  new HashMap<String,IntervalTree<String>>();
							String distantRegionName = splitStrings[0];
							String distantRegion = splitStrings[1];
							distantLocString.add(distantRegionName);
							
							BufferedReader br = new BufferedReader(new FileReader(distantRegion));
							String line;
							while( (line = br.readLine()) != null){
								if(line.startsWith("#"))
									continue;
								String[] splitLines = line.split("\t");
								if(splitLines.length<3){
									continue;
								}else{
									String chr = splitLines[0];
									int start = Integer.parseInt(splitLines[1]);
									int end = Integer.parseInt(splitLines[2]);
									IntervalTree<String> tree;
									if(distantLocCollections.containsKey(chr)){
										tree = distantLocCollections.get(chr);
									}else{
										tree = new IntervalTree<String>();
									}
									String strand = ".";
									if(splitLines.length >= 6){
										if(splitLines[5].equalsIgnoreCase("-")){
											strand = "-";
										}else if(splitLines[5].equalsIgnoreCase("+")){
											strand = "+";
										}
										//strand = splitLines[5].equalsIgnoreCase("-") ? "-" : "+";
									}
									tree.put(start, end, strand);
									distantLocCollections.put(chr, tree);
								}
								
								
								
								
							}
							br.close();
							distantLocStringCollections.put(distantRegionName, distantLocCollections);
						
						}
					}
					
					
					HashMap<String, Pair<Integer, TabixFeatureReader<BEDFeature, ?>>> valueBedReaders = null;
					LinkedHashSet<String> valueBedLocString = new LinkedHashSet<String>();
					if(valueBeds != null){
						log.info("Loading value interval bed file ... ");
						valueBedReaders =  new HashMap<String, Pair<Integer, TabixFeatureReader<BEDFeature, ?>>>();
						for(String valueBedString : valueBeds){
							String[] splitStrings = valueBedString.split(":");
							if(splitStrings.length < 3){
								throw new IllegalArgumentException("need to provide trackname:extRegion:trackFile for valueBeds");
							}
							String valueBedName = splitStrings[0];
							int valueBedExt = Integer.parseInt(splitStrings[1]);
							String valueRegion = splitStrings[2];
							valueBedLocString.add(valueBedName);
							valueBedReaders.put(valueBedName,new Pair<Integer, TabixFeatureReader<BEDFeature, ?>>(valueBedExt,new TabixFeatureReader(valueRegion, new BEDCodec())));
						}
						
					}
					
					HashMap<String, Pair<Integer, BigWigFileReader>> valueWigReaders = null;
					LinkedHashSet<String> valueWigLocString = new LinkedHashSet<String>();
					if(valueWigs != null){
						log.info("Loading value interval big wig file ... ");
						valueWigReaders =  new HashMap<String, Pair<Integer,BigWigFileReader>>();
						for(String valueWigString : valueWigs){
							String[] splitStrings = valueWigString.split(":");
							if(splitStrings.length < 3){
								throw new IllegalArgumentException("need to provide trackname:extRegion:trackFile for valueWigs");
							}
							String valueWigName = splitStrings[0];
							int valueWigExt = Integer.parseInt(splitStrings[1]);
							String valueRegion = splitStrings[2];
							valueWigLocString.add(valueWigName);
							valueWigReaders.put(valueWigName,new Pair<Integer, BigWigFileReader>(valueWigExt,new BigWigFileReader(new File(valueRegion).toPath())));
						}
						
					}
					
					LinkedHashSet<String> kmerCollections = new LinkedHashSet<String>();
					if(kmerString != null){
						log.info("Loading selected K-mer file ... ");	
						BufferedReader br = new BufferedReader(new FileReader(kmerString));
						String line;
						while( (line = br.readLine()) != null){
							if(line.startsWith("#"))
								continue;
							kmerCollections.add(line);
							
							
						}
						br.close();
					}else{
						log.info("Automate generate all k-mer until length " + kmerLen);
						for(int i = 2; i <=  kmerLen; i++){
							for(byte[] kmer : SequenceUtil.generateAllKmers(i)){
								kmerCollections.add(new String(kmer));
								
							}
						}
						
					}
					
					String header = "";
					if(overlapLocString.size()>0){
						for(String key : overlapLocString){
							header = header + "\t" + key;
						}
					}
					if(distantLocString.size()>0){
						for(String key : distantLocString){
							header = header + "\t" + key;
						}
					}
					if(valueBedLocString.size()>0){
						for(String key : valueBedLocString){
							header = header + "\t" + key;
						}
					}
					if(valueWigLocString.size()>0){
						for(String key : valueWigLocString){
							header = header + "\t" + key;
						}
					}
					if(kmerCollections.size()>0){
						//if(useRefSeqBaseKmer){
							for(String key : kmerCollections){
								header = header + "\t" + key;
							}
						//}else{
							
						//}
						
					}
					
					
					log.info("Loading CpG interval file ... ");
					HashMap<String,IntervalTree<String>> cpgCollections = new HashMap<String,IntervalTree<String>>();
						GZIPInputStream gzipInputStream1 = null;
						BufferedReader br;
						if(allCpgFile.endsWith(".gz")){
							gzipInputStream1 = new GZIPInputStream(new FileInputStream(cpgListFile));
							br = new BufferedReader(new InputStreamReader(gzipInputStream1));
							
						}else{
							br = new BufferedReader(new FileReader(cpgListFile));
						}
							
							String line;
							
							while( (line = br.readLine()) != null){
								if(line.startsWith("#"))
									continue;
								String[] splitLines = line.split("\t");
								if(splitLines.length<3){
									continue;
								}
								String chr = splitLines[0];
								int start = Integer.parseInt(splitLines[1]);
								int end = Integer.parseInt(splitLines[2]);
								IntervalTree<String> tree;
								
								if(cpgCollections.containsKey(chr)){
									tree = cpgCollections.get(chr);
								}else{
									tree = new IntervalTree<String>();
								}
								String strand = ".";
								if(splitLines.length >= 6){
									if(splitLines[5].equalsIgnoreCase("-")){
										strand = "-";
									}else if(splitLines[5].equalsIgnoreCase("+")){
										strand = "+";
									}
									//strand = splitLines[5].equalsIgnoreCase("-") ? "-" : "+";
								}
								tree.put(start, end, strand);
								cpgCollections.put(chr, tree);
							}
							if(cpgListFile.endsWith(".gz")){
								gzipInputStream1.close();
							}
							br.close();

					log.info("Output value for each CpG  ... ");
					FileOutputStream output = new FileOutputStream(detailFile);
					OutputStreamWriter writer = new OutputStreamWriter(new GZIPOutputStream(output), "UTF-8");
					
					//basic
					writer.write("chr\tstart\tend");
					if(!skipCpgDist){
						writer.write("\tdist_nearest_CpG");
					}
					writer.write(header + "\n");
					//fragment k-mer
					
					
					//
					
					long i = 0;
					String prevChr = "";
					for(String chr : cpgCollections.keySet()){
						if(chr.equalsIgnoreCase("chrM")){
							continue;
							
						}
						if(prevChr.equalsIgnoreCase("")){
							prevChr = chr;
							refParser.setCurrentSequence(chr);
						}
						
						IntervalTree<String> cpgChrCollections = cpgCollections.get(chr);
						Iterator<Node<String>> cpgIterator = cpgChrCollections.iterator();
						while(cpgIterator.hasNext()){
							Node<String> cpg = cpgIterator.next();
							int start = cpg.getStart();
							int end = cpg.getEnd();
							int fragMostLeft = start+1;
							int fragMostRight = end;
							
							
						if(!chr.equalsIgnoreCase(prevChr)){
							refParser.close();
							refParser.setCurrentSequence(chr);
							prevChr = chr;
						}
						byte[] refBasesExt = CcInferenceUtils.toUpperCase(refParser.loadFragment(end-1-kmerExt, kmerExt*2+1).getBytes());
						byte refBase = refBasesExt[kmerExt];
						
						//if(end==16217220){
						//	log.info(new String(refBasesExt) + "\t" + start + "\t" + end + "\t" + (end-1-kmerExt) + "\t" + (kmerExt*2+1) + "\t" + kmerExt + "\t" + (char)refBasesExt[kmerExt]);
						//	
						//}
						
						
						HashMap<String, Double> kmerMapsRef = new HashMap<String, Double>();
						//if(useRefSeqBaseKmer){
							for(int j = 2; j <= kmerLen; j++){
								kmerMapsRef.putAll(CcInferenceUtils.kmerFreqSearch(refBasesExt, j));
								
							}
						//}
						
						//nearest cpg's distance in reference genome
							
						IntervalTree<String> cpgLocCollections = null;
						double nearestCpg = Double.NaN;
						if(!skipCpgDist){
							cpgLocCollections = allCpgLocCollections.get(chr);
							Iterator<Node<String>> upstreamCpgIt = null;
							Iterator<Node<String>> downstreamCpgIt = null;
							if(BaseUtils.basesAreEqual(refBase, BaseUtilsMore.C)){
								upstreamCpgIt = cpgLocCollections.reverseIterator(start-1, end-1);
								downstreamCpgIt = cpgLocCollections.iterator(start+2, end+2);
								
							}else if(BaseUtils.basesAreEqual(refBase, BaseUtilsMore.G)){
								upstreamCpgIt = cpgLocCollections.reverseIterator(start-2, end-2);
								downstreamCpgIt = cpgLocCollections.iterator(start+1, end+1);
							}else{
								continue;
							}
							
							if(upstreamCpgIt== null || !upstreamCpgIt.hasNext()){
								IntervalTree.Node<String> downstream = downstreamCpgIt.next();
								nearestCpg = CcInferenceUtils.intervalDistance(downstream,cpg);
								//System.err.println(downstream.toString());
							}else if(downstreamCpgIt==null || !downstreamCpgIt.hasNext()){
								IntervalTree.Node<String> upstream = upstreamCpgIt.next();
								nearestCpg = CcInferenceUtils.intervalDistance(upstream,cpg);
								//System.err.println(upstream.toString());
							}else{
								IntervalTree.Node<String> upstream = upstreamCpgIt.next();
								IntervalTree.Node<String> downstream = downstreamCpgIt.next();
								//System.err.println(upstream.toString());
								//System.err.println(downstream.toString());
								
								int dist1 = CcInferenceUtils.intervalDistance(upstream, cpg);
								int dist2 = CcInferenceUtils.intervalDistance(downstream, cpg);
								if(Math.abs(dist1) < Math.abs(dist2)){
									nearestCpg = dist1;
								}else{
									nearestCpg = dist2;
								}
							}
						}
						
						//overlap with feature in reference genome
						
						HashMap<String, Integer> overlapStatCollections = new HashMap<String, Integer>();
						if(overlapRegions!=null && !overlapRegions.isEmpty()){
							for(String key : overlapLocStringCollections.keySet()){
								HashMap<String,IntervalTree<Integer>> tmp = overlapLocStringCollections.get(key);
								if(tmp.containsKey(chr)){
									if(tmp.get(chr).minOverlapper(start, end)==null){
										overlapStatCollections.put(key, 0);
									}else{
										overlapStatCollections.put(key, 1);
									}
									
								}else{
									overlapStatCollections.put(key,0);
								}
								
							}
						}
						
						//distance with feature in reference genome
						HashMap<String, Integer> distStatCollections = new HashMap<String, Integer>();
						if(distantRegions!=null && !distantRegions.isEmpty()){
							for(String key : distantLocStringCollections.keySet()){
								IntervalTree<String> locCollections = distantLocStringCollections.get(key).get(chr);
								
								//IntervalTree.Node<String> upstream = locCollections.max(start, end);
								//IntervalTree.Node<String> downstream = locCollections.min(start, end);
								int distanceNearest = Integer.MAX_VALUE;
								if(locCollections!=null && locCollections.size()>0){
									Iterator<Node<String>> upstreamIt = locCollections.reverseIterator(start, end);
									Iterator<Node<String>> downstreamIt = locCollections.iterator(start, end);
									if(!upstreamIt.hasNext()){
										IntervalTree.Node<String> downstream = locCollections.min(start, end);
										distanceNearest = CcInferenceUtils.intervalDistance(downstream,cpg);
										//System.err.println(downstream.toString());
									}else if(!downstreamIt.hasNext()){
										IntervalTree.Node<String> upstream = locCollections.max(start, end);
										distanceNearest = CcInferenceUtils.intervalDistance(upstream,cpg);
										//System.err.println(upstream.toString());
									}else{
										IntervalTree.Node<String> upstream = locCollections.max(start, end);
										IntervalTree.Node<String> downstream = locCollections.min(start, end);
										//System.err.println(upstream.toString());
										//System.err.println(downstream.toString());
										
										int dist1 = CcInferenceUtils.intervalDistance(upstream, cpg);
										int dist2 = CcInferenceUtils.intervalDistance(downstream, cpg);
										if(Math.abs(dist1) < Math.abs(dist2)){
											distanceNearest = dist1;
										}else{
											distanceNearest = dist2;
										}
									}
								}
								
								
								distStatCollections.put(key, distanceNearest);
							}
						}
						
						//value in bed file
						HashMap<String, Double> valBedStatCollections = new HashMap<String, Double>();
						if(valueBedReaders != null){
							for(String key : valueBedReaders.keySet()){
								int range = valueBedReaders.get(key).getFirst();
								boolean mean0 = range < 0 ? true : false;
								if(mean0){
									range = 0-range;
								}
								TabixFeatureReader<BEDFeature, ?> bedReader = valueBedReaders.get(key).getSecond();
								CloseableTribbleIterator<BEDFeature> featureIt = bedReader.query(chr, (start-range < 0 ? 1 : start-range+1), end+range);
								DescriptiveStatistics statFeature = new DescriptiveStatistics();
								while(featureIt.hasNext()){
									BEDFeature term = featureIt.next();
									if(!Double.isNaN(term.getScore())){
										statFeature.addValue(term.getScore());
									}
									
								}
								featureIt.close();
								if(statFeature.getN()>0){
									if(mean0){
										valBedStatCollections.put(key, (double)statFeature.getSum()/(double)(range*2+1));
									}else{
										valBedStatCollections.put(key, statFeature.getMean());
									}
									
								}else{
									valBedStatCollections.put(key, Double.NaN);
								}
								
							}
						}
						
						
						//value in wig file
						HashMap<String, Double> valWigStatCollections = new HashMap<String, Double>();
						if(valueWigReaders != null){
							for(String key : valueWigReaders.keySet()){
								int range = valueWigReaders.get(key).getFirst();
								if(range < 0){
									BigWigFileReader wigReader = valueWigReaders.get(key).getSecond();
									range = 0-range;
									SummaryStatistics statFeature = wigReader.queryStats(chr, (start-range < 0 ? 1 : start-range), end+range);
									
									if(statFeature.getN()>0){
										valWigStatCollections.put(key, (double)statFeature.getSum()/(double)(range*2+1));
									}else{
										valWigStatCollections.put(key, Double.NaN);
									}
									
								}else{
									BigWigFileReader wigReader = valueWigReaders.get(key).getSecond();
									SummaryStatistics statFeature = wigReader.queryStats(chr, start-range, end+range);
									if(statFeature.getN()>0){
										valWigStatCollections.put(key, statFeature.getMean());
									}else{
										valWigStatCollections.put(key, Double.NaN);
									}
								}
								
								
								
							}
						}
						
						
							
							
								//System.err.println(CcInferenceUtils.getFragOffsetFromReadsOffset(r, offSet));
							writer.write(chr + "\t" + start + "\t" + end);
							if(!skipCpgDist){
								writer.write("\t" + nearestCpg);
							}
							//overlap regions
							if(overlapStatCollections.size()>0){
								for(String key : overlapLocString){
									writer.write("\t" + overlapStatCollections.get(key));
								}
							}
							
							//distant regions
							if(distStatCollections.size() > 0){
								for(String key : distantLocString){
									writer.write("\t" + distStatCollections.get(key));
								}
							}
							
							//valBed regions
							if(valBedStatCollections.size()>0){
								for(String key : valueBedLocString){
									writer.write("\t" + String.format("%.3f",valBedStatCollections.get(key)));
								}
							}
							
							//valWig regions
							if(valWigStatCollections.size()>0){
								for(String key : valueWigLocString){
									writer.write("\t" + String.format("%.3f",valWigStatCollections.get(key)));
								}
							}
							//k-mer in reference genome
							//if(kmerMapsRef.size()>0 && useRefSeqBaseKmer){
							if(kmerMapsRef.size()>0 ){
								for(String key : kmerCollections){
									writer.write("\t" + String.format("%.3f",kmerMapsRef.get(key)));
								}
							}
							
							
							
							writer.write("\n");
							points++;
							
							
						}
						i++;
						if(i % 1000 == 0){
							log.info("Processing Cpg " + i + " ...");
							writer.flush();
						}
						
					}
					writer.close();
					output.close();
					
					

					refParser.closeParser();;

					
					if(valueBedReaders != null){
						for(String key : valueBedReaders.keySet()){
							valueBedReaders.get(key).getSecond().close();
						}
					}
					if(valueWigReaders != null){
						for(String key : valueWigReaders.keySet()){
							valueWigReaders.get(key).getSecond().close();
						}
					}
					
					finish();

	}
	
	

	

	
	
	private void initiate(){
		startTime = System.currentTimeMillis();
	}

	private void finish(){
		long endTime   = System.currentTimeMillis();
		double totalTime = endTime - startTime;
		totalTime /= 1000;
		double totalTimeMins = totalTime/60;
		double totalTimeHours = totalTime/3600;
		
		log.info("Counted " + points + " data points in total");
		log.info("CpgMultiMetricsStatsNoBam's running time is: " + String.format("%.2f",totalTime) + " secs, " + String.format("%.2f",totalTimeMins) +  " mins, " + String.format("%.2f",totalTimeHours) +  " hours");
	}
	

}
