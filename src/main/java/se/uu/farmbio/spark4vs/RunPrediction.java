/*** RunPrediction.java ***/

package se.uu.farmbio.spark4vs;

import java.io.IOException;
import java.io.StringReader;
import java.util.List;
import java.util.Arrays;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import se.uu.farmbio.spark4vs.model.SignLibsvmModel;
import se.uu.farmbio.spark4vs.model.SignLibsvmUtils;
import se.uu.farmbio.spark4vs.io.SDFInputFormat;

import scala.Tuple2;
import spark.api.java.*;
import spark.api.java.function.*;
import spark.storage.StorageLevel;

import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.conf.Configuration;

/**
 * A class to filter out molecules by predictive models.
 * 
 * @author ola and laeeq
 *
 */
 
public class RunPrediction{
	
	public static void main(String[] args) throws IOException{
	
		if (args.length < 2) {
      		System.err.println("Usage: SparkForVirtualScreening <master> <SDFfile>");
      		System.exit(1);
   	 	}
   	 	
   	 	//Setting system properties
   	 	System.setProperty("spark.cores.max", "7");
		System.setProperty("spark.default.parallelism", "40");
		//System.setProperty("spark.rdd.compress","true");
		
		JavaSparkContext jsc = new JavaSparkContext(args[0], "Spark for Virtual Screening", 
		"$YOUR_SPARK_HOME", "target/spark4vs-1.0-jar-with-dependencies.jar");
		
		
		Configuration conf = new Configuration();
		
		//Basic RDD with the complete file in it
		JavaPairRDD<LongWritable, Text> molecules = jsc.newAPIHadoopFile(args[1], SDFInputFormat.class, LongWritable.class, Text.class, conf);
		
		molecules.cache();
		
		String amesModelFile = RunPrediction.class.getResource("/ames.model").getFile();
		String amesSignFile = RunPrediction.class.getResource("/ames.sign").getFile();
		String hergModelFile = RunPrediction.class.getResource("/hERG@PKKB-all-data-classi.model").getFile();
		String hergSignFile = RunPrediction.class.getResource("/hERG@PKKB-all-data-classi.sign").getFile();

		//Initialize the models into memory
		final SignLibsvmModel amesModel = SignLibsvmUtils.ModelFromFile(amesModelFile,amesSignFile);
		final SignLibsvmModel hergModel = SignLibsvmUtils.ModelFromFile(hergModelFile,hergSignFile);
		
		//Each Molecule predicted for mutagenicity and inhibition and returns value(1 or 0) for all molecules
		JavaRDD<Integer> winners = molecules.map(new Function<Tuple2<LongWritable,Text>,Integer>(){
				@Override
				public Integer call(Tuple2<LongWritable,Text> pair) throws CDKException, IOException
				{
					
					//Reading molecule from RDD's Second Argument(value in (key,value) pair)
					StringReader sreader = new StringReader(pair._2().toString());
					MDLV2000Reader reader = new MDLV2000Reader(sreader);
					ChemFile chemFile = (ChemFile) reader.read((ChemObject) new ChemFile());
                	List<IAtomContainer> containersList = ChemFileManipulator.getAllAtomContainers(chemFile);
                	IAtomContainer mol = containersList.get(0);
                
					//Calculating signature for molecules
					List<String> signatures=SignLibsvmUtils.calculateSignatures(mol, 0, 3);
					
					//Predicting winner molecules by comparing with svm model
					double amesRes = amesModel.predict(signatures);
					double hergRes = hergModel.predict(signatures);
		
					//Predict mutagenicity = 0 and Predict herg inhibition = 1
					if ((amesRes!=0) && (hergRes!=1))
						return 1;
					else
						return 0;	
				}
		});
    
    winners.cache();
    
    //Performing Sum of values to find winner molecules
    Integer total = winners.reduce(new Function2<Integer, Integer, Integer>() {
      @Override
      public Integer call(Integer a, Integer b) {
        return a + b;
      }
    });
  	
	System.out.println("Predicted " + total + " mols in " + winners.count());  
	System.exit(0);		
	}
		
}
