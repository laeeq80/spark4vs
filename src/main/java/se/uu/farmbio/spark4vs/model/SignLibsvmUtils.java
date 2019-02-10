package se.uu.farmbio.spark4vs.model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import libsvm.svm;
import libsvm.svm_model;

import org.apache.log4j.Logger;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.signature.AtomSignature;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;


public class SignLibsvmUtils {

	private static final Logger logger = Logger.getLogger(SignLibsvmUtils.class);

	public static SignLibsvmModel ModelFromFile(String modelFile, String signaturesFile)
			throws IOException{

		List<String> signatures=readSignaturesFile(signaturesFile);
		svm_model svmModel = svm.svm_load_model(modelFile);
		
		SignLibsvmModel model = new SignLibsvmModel(svmModel,signatures);

		return model;
	}




	/**
	 * Read a list of signatures from file into an arraylist
	 * 
	 * @param signaturesPath
	 * @return
	 * @throws IOException 
	 * @throws DSException
	 */
	public static List<String> readSignaturesFile(String signaturesPath) throws IOException {

		logger.debug("Reading signature file: " + signaturesPath);

		List<String> signatures = new ArrayList<String>(); // Contains modelSignatures. We use the indexOf to retrieve the order of specific modelSignatures in descriptor array.
		BufferedReader signaturesReader=null;
		try {
			signaturesReader = new BufferedReader(new FileReader(new File(signaturesPath)));
			String signature;
			while ( (signature = signaturesReader.readLine()) != null ) {
				//				if (modelSignatures.contains(signature))
				//					throw new DSException("Duplicate signature in modelSignatures list");
				signatures.add(signature);
			}
		}finally{
				if (signaturesReader!=null)
					signaturesReader.close();
		}

		logger.debug("Reading signature file: " + signaturesPath + " completed successfully");

		return signatures;

	}

	public static List<String> calculateSignatures(IAtomContainer mol, int startheight, int endheight) throws CDKException{

		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		CDKHueckelAromaticityDetector.detectAromaticity(mol);
		mol=AtomContainerManipulator.removeHydrogens(mol);
		
		List<String> ret = new ArrayList<String>();

		for (int h=startheight; h< endheight; h++){

			for (int i=0; i<mol.getAtomCount(); i++){
				AtomSignature as = new AtomSignature(i, h, mol);
				ret.add(as.toCanonicalString());
			}
		}

		return ret;
		
	}

}
