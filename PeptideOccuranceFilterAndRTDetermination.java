/**
 * Filter peptide by number of appearance, get average RT
 * @author Bruce Lab
 *
 */
import java.io.*;
import java.util.*;

public class PeptideOccuranceFilterAndRTDetermination {
	public static void main(String[] args) {
	}
	public static void filter_largeRTRange(String in, String out) {
		try {
			FileReader fr = new FileReader(in);
			BufferedReader br = new BufferedReader(fr);
			FileOutputStream fout = new FileOutputStream(out);
			PrintStream ps = new PrintStream(fout);
			String line = br.readLine();
			while (line != null) {
				String[] arr = line.split("\t");
				String[] rts = arr[3].split("-");
				double diff = Double.valueOf(rts[1]) - Double.valueOf(rts[0]);
				if (diff < 2) {
					double rt = (Double.valueOf(rts[1]) + Double.valueOf(rts[0])) / 2;
					ps.println(arr[0] + "\t" + arr[1] + "\t" + arr[2] + "\t" + rt + "\t" + arr[4]);
				}
				line = br.readLine();
			}
			ps.close();
			fout.close();
			br.close();
			fr.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
	}
	public static void filter_new(ArrayList<Peptide>peps, int minOccur, String out) {
		try {
			FileOutputStream fout = new FileOutputStream(out);
			PrintStream ps = new PrintStream(fout);
			for (Peptide p : peps) {
				if (p.files.size() < minOccur) continue;
				String rt = rtRange(p.rt);
				ps.println(p.pname + "\t" + p.seq + "\t" + p.mz + "\t" + rt + "\t" + p.multiChargeFlag);
			}
			ps.close();
			fout.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
	}
	public static void filter(ArrayList<Peptide>peps, int minOccur, String out) {
		try {
			FileOutputStream fout = new FileOutputStream(out);
			PrintStream ps = new PrintStream(fout);
			for (Peptide p : peps) {
				if (p.files.size() < minOccur) continue;
				double rt = getRT(p.rt);
				ps.println(p.pname + "\t" + p.seq + "\t" + p.mz + "\t" + rt + "\t" + p.multiChargeFlag);
			}
			ps.close();
			fout.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
	}
	public static String rtRange(ArrayList<Double> rt) {
		String rtRange = "";
		double min = rt.get(0);
		double max = 0;
		for (double d : rt) {
			if (d < min) min = d;
			if (d > max) max = d;
		}
		rtRange = min + "-" + max;
		return rtRange;
	}
	public static double getRT(ArrayList<Double> rt) {
		double sum = 0; 
		for (double d : rt) sum += d;
		return sum / rt.size();
	}
	public static ArrayList<Peptide> assignMultiChargeFlag(ArrayList<Peptide> peps) {
		Set<Integer> visited = new HashSet<Integer>();
		for (int i = 0; i < peps.size() - 1; i++) {
			if (visited.contains(i)) continue;
			ArrayList<Integer> index = new ArrayList<Integer>();
			index.add(i);
			String seq = peps.get(i).seq;
			int j = i + 1;
			while (j < peps.size()) {
				String seqtmp = peps.get(j).seq;
				if (seq.equals(seqtmp)) index.add(j);
				j++;
			}
			if (index.size() > 1) {
				for (int k : index) {
					peps.get(k).multiChargeFlag = true;
					visited.add(k);
				}
			} else {
				visited.add(i);
			}
		}
		return peps;
	}
	public static ArrayList<Peptide> getPeptide(String in) {
		ArrayList<Peptide> peptides = new ArrayList<Peptide>();
		Map<String, Peptide> pepMap = new HashMap<String, Peptide>();
		try {
			FileReader fr = new FileReader(in);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			while (line != null) {
				String[] arr = line.split("\t");
				String key = arr[2] + arr[3];
				if (pepMap.containsKey(key)) {
					Peptide pep = pepMap.get(key);
					pep.files.add(arr[0]);
					pep.rt.add(Double.valueOf(arr[4]));
				} else {
					Peptide pep = new Peptide();
					pep.pname = arr[1];
					pep.seq = arr[2];
					pep.mz = Double.valueOf(arr[3]);
					pep.files.add(arr[0]);
					pep.rt.add(Double.valueOf(arr[4]));
					pepMap.put(key, pep);
					peptides.add(pep);
				}
				line = br.readLine();
			}
			br.close();
			fr.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		return peptides;
	}
}
class Peptide {
	public String pname;
	public String seq;
	public int charge;
	public double mz;
	public ArrayList<Double> rt;
	public ArrayList<String> files;
	public boolean multiChargeFlag;
	public Peptide() {
		rt = new ArrayList<Double>();
		files = new ArrayList<String>();
		multiChargeFlag = false;
	}
}
