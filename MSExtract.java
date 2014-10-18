/**
 * Mass spec peak extractor
 * @author brucelab
 *
 */
import java.io.*;
import org.systemsbiology.jrap.stax.*;
import java.util.*;

public class MSExtract {
	public static void main(String[] args) {
	}
	@SuppressWarnings("unchecked")
	public static void getPeakAreaWrapper_fast(String filelist, String in, double rtTolerance, double runLen, double ppmTolerance) {
		//////////////////////////////////////////Read in all the precursors
		ArrayList<Precursor> precursors = new ArrayList<Precursor>();
		try {
			FileReader fr = new FileReader(in);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			while (line != null) {
				String[] arr = line.split("\t");
				Precursor p = new Precursor();
				p.mz = Double.valueOf(arr[2]);
				p.lb = Double.valueOf(arr[3]) - rtTolerance;
				p.ub = Double.valueOf(arr[3]) + rtTolerance;
				precursors.add(p);
				line = br.readLine();
			}
			br.close();
			fr.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		if (precursors.isEmpty()) return;
		////////////////////////////////////////////Sort precursor, generate empty array to store areas
		Collections.sort(precursors);
		int numPre = precursors.size();
		double[] areas = new double[numPre];
		///////////////////////////////////////////read in file path and name
		ArrayList<String> files = new ArrayList<String>();
		try {
			FileReader fr = new FileReader(filelist);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			while (line != null) {
				files.add(line);
				line = br.readLine();
			}
			br.close();
			fr.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		/////////////////////////////////////////extract peak area
		try {
			FileOutputStream fout = new FileOutputStream(in + ".out");
			PrintStream ps = new PrintStream(fout);
			for (String f : files) {
				System.out.println(f);
				MSXMLParser parser = new MSXMLParser(f);
				for (int i = 1; i <= parser.getMaxScanNumber(); i++) {
					System.out.println(i);
					if (parser.rap(i).header.getMsLevel() != 1) continue;
					Scan scan = parser.rap(i);
					double rt = parseRT(parser.rapHeader(i).getRetentionTime());
					double[][] peaklist = scan.getMassIntensityList();
					for (int j  = 0; j < precursors.size(); j++) {
						Precursor  p = precursors.get(j);
						System.out.println(p.mz);
						if (rt < p.lb || rt > p.ub) continue;
						areas[j] += getIntensity(peaklist, p.mz, ppmTolerance);
					}
				}
				for (int i = 0; i < precursors.size(); i++) {
					ps.println(f + "\t" + precursors.get(i).mz + "\t" + areas[i]);
					areas[i] = 0;
				}
			}
			ps.close();
			fout.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
	}
	public static double getIntensity(double[][] peaklist, double mz, double ppmtolerance) {
		double intensity = 0;
		double lb = mz - ppmtolerance / 1000000 * mz;
		double ub = mz + ppmtolerance / 1000000 * mz;
		for (int i = 0; i < peaklist[0].length; i++) {
			if (peaklist[0][i] >= lb && peaklist[0][i] <= ub) {
				intensity += peaklist[1][i];
			}
		}
		return intensity;
	}
	public static void getPeakAreaWrapper_new(String filelist, String in, double rtTolerance) {
		ArrayList<String> flist = new ArrayList<String>();
		try {
			FileReader fr = new FileReader(filelist);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			while (line != null) {
				flist.add(line);
				line = br.readLine();
			}
			br.close();
			fr.close();
		} catch (IOException e) {
			System.err.println(e.getMessage());
		}
		try {
			FileReader fr = new FileReader(in);
			BufferedReader br = new BufferedReader(fr);
			FileOutputStream fout = new FileOutputStream(in + ".out");
			PrintStream ps = new PrintStream(fout);
			String line = br.readLine();
			while (line != null) {
				String[] arr = line.split("\t");
				for (String s : flist) {
					String fileName = s;
					double rt = Double.valueOf(arr[3]);
					double mz = Double.valueOf(arr[2]);
					System.out.println(fileName + "\t" + rt + "\t" + mz);                             /////test
					double area = 0;
					double rtlb = rt - rtTolerance;
					double rtub = rt + rtTolerance;
					MSXMLParser parser = new MSXMLParser(fileName);
					int mxScan = parser.getMaxScanNumber();
					int currentScan = (int) (rtlb / 160 * mxScan);
					ScanHeader currentScanHeader = parser.rapHeader(currentScan);
					double currentRT = parseRT(currentScanHeader.getRetentionTime());
					while (currentRT > rtlb) {
						currentScan -= 100;
						currentScanHeader = parser.rapHeader(currentScan);
						currentRT = parseRT(currentScanHeader.getRetentionTime());
					}
					while (currentRT < rtub) {
						//System.out.println(currentScan);                                             /////test
						if (currentRT < rtlb) {
							currentScan++;
							currentScanHeader = parser.rapHeader(currentScan);
							currentRT = parseRT(currentScanHeader.getRetentionTime());
							continue;
						}
						if (currentScanHeader.getMsLevel() != 1) {
							currentScan++;
							currentScanHeader = parser.rapHeader(currentScan);
							currentRT = parseRT(currentScanHeader.getRetentionTime());
							continue;
						}
						Scan scan = parser.rap(currentScan);
						double[][] list = scan.getMassIntensityList();
						for (int i = 0; i < list[0].length; i++) {
							if (list[0][i] > mz + mz * 10 / 1000000) break;
							if (list[0][i] < mz - mz * 10 / 1000000) continue;
							area += list[1][i];
						}
						currentScan++;
						currentScanHeader = parser.rapHeader(currentScan);
						currentRT = parseRT(currentScanHeader.getRetentionTime());
					}					
					ps.println(s + "\t" + line + "\t" + area);
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
	public static double parseRT(String in) {
		return Double.valueOf(in.substring(2, in.length() - 1)) / 60;
	}
	public static void getPeakAreaWrapper(String in, double rtTolerance) {
		try {
			FileReader fr = new FileReader(in);
			BufferedReader br = new BufferedReader(fr);
			FileOutputStream fout = new FileOutputStream(in + ".out");
			PrintStream ps = new PrintStream(fout);
			String line = br.readLine();
			while (line != null) {
				String[] arr = line.split("\t");
				String fileName = arr[0];
				double rt = Double.valueOf(arr[6]);
				double mz = Double.valueOf(arr[3]);
				System.out.println(fileName + "\t" + rt + "\t" + mz);                             /////test
				double area = 0;
				double rtlb = rt - rtTolerance;
				double rtub = rt + rtTolerance;
				MSXMLParser parser = new MSXMLParser(fileName);
				int mxScan = parser.getMaxScanNumber();
				int currentScan = (int) (rtlb / 160 * mxScan);
				ScanHeader currentScanHeader = parser.rapHeader(currentScan);
				double currentRT = parseRT(currentScanHeader.getRetentionTime());
				while (currentRT > rtlb) {
					currentScan -= 100;
					currentScanHeader = parser.rapHeader(currentScan);
					currentRT = parseRT(currentScanHeader.getRetentionTime());
				}
				while (currentRT < rtub) {
					System.out.println(currentScan);                                             /////test
					if (currentRT < rtlb) {
						currentScan++;
						currentScanHeader = parser.rapHeader(currentScan);
						currentRT = parseRT(currentScanHeader.getRetentionTime());
						continue;
					}
					if (currentScanHeader.getMsLevel() != 1) {
						currentScan++;
						currentScanHeader = parser.rapHeader(currentScan);
						currentRT = parseRT(currentScanHeader.getRetentionTime());
						continue;
					}
					Scan scan = parser.rap(currentScan);
					double[][] list = scan.getMassIntensityList();
					for (int i = 0; i < list[0].length; i++) {
						if (list[0][i] > mz + mz * 10 / 1000000) break;
						if (list[0][i] < mz - mz * 10 / 1000000) continue;
						area += list[1][i];
					}
					currentScan++;
					currentScanHeader = parser.rapHeader(currentScan);
					currentRT = parseRT(currentScanHeader.getRetentionTime());
				}
				
				ps.println(line + "\t" + area);
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
	public static ArrayList<Integer> getScanListFromRT(String fname, double rt, double rtTolerance) {
		ArrayList<Integer> tmpScanList = new ArrayList<Integer>();
		MSXMLParser parser = new MSXMLParser(fname);
		double rtlb = rt - rtTolerance;
		double rtub = rt + rtTolerance;
		int scanNum = 1;
		double tmprt = getRT(fname, scanNum);
		while (tmprt < rtub) {
			if (tmprt > rtlb)  {
				if (parser.rap(scanNum).header.getMsLevel() == 1) tmpScanList.add(scanNum);
			}
			//System.out.println(scanNum);
			scanNum++;
			tmprt = getRT(fname, scanNum);
		}
		return tmpScanList;		
	}
	public static void getRTwrapper(String in) {
		try {
			FileReader fr = new FileReader(in);
			BufferedReader br = new BufferedReader(fr);
			FileOutputStream fout = new FileOutputStream(in + ".out");
			PrintStream ps = new PrintStream(fout);
			String line = br.readLine();
			while (line != null) {
				String[] arr = line.split("\t");
				ps.println(line + "\t" + getRT(arr[0], Integer.valueOf(arr[4].trim())));
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
	public static double getRT(String fname, int scanNum) {
		double rt = 0;
		MSXMLParser msxmlParser = new MSXMLParser(fname);
		Scan scan = msxmlParser.rap(scanNum);
		ScanHeader scanHeader = scan.getHeader();
		String rtTxt = scanHeader.getRetentionTime();
		rt = Double.valueOf(rtTxt.substring(2, rtTxt.length() - 1)) / 60;
		return rt;
	}
	public static void mxScanWrapper(String in) {
		try {
			FileReader fr = new FileReader(in);
			BufferedReader br = new BufferedReader(fr);
			FileOutputStream fout = new FileOutputStream(in + ".out");
			PrintStream ps = new PrintStream(fout);
			String line = br.readLine();
			while (line != null) {
				ps.println(line + "\t" + returnMxScan(line));
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
	public static int returnMxScan(String in) {
		MSXMLParser fileParser = new MSXMLParser(in);
		return fileParser.getMaxScanNumber();
	}
	public static void wrapper(String in, String out, String mzXMLdir) {
		try {
			FileReader fr = new FileReader(in);
			BufferedReader br = new BufferedReader(fr);
			FileOutputStream fout = new FileOutputStream(out);
			PrintStream ps = new PrintStream(fout);
			String line = br.readLine();
			while (line != null) {
				String[] arr = line.split("\t");
				String mzXMLFileName = mzXMLdir + "/" + arr[0] + ".mzXML";
				double mz = Double.valueOf(arr[3]);
				int scanNum = Integer.valueOf(arr[4]);
				int[] scanlist = getScanList(scanNum, mzXMLFileName);
				double area = integratePeakBoxCar(scanlist, mz, mzXMLFileName);
				ps.println(line + "\t" + area);
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
	public static double integratePeak(int[] scanlist, double mz, String fileName) {
		double peakArea = 0;
		double lB, uB;
		lB = mz - 0.000020 * mz;
		uB = mz + 0.000020 * mz;
		MSXMLParser parser = new MSXMLParser(fileName);
		for (int scanNum : scanlist) {
			Scan scan = parser.rap(scanNum);
			double[][] massIntensityList = scan.getMassIntensityList();
			for(int i = 0; i < massIntensityList[0].length; i++) {
				if (massIntensityList[0][i] < lB) continue;
				if (massIntensityList[0][i] > uB) break;
				peakArea += massIntensityList[1][i];
			}
		}
		return peakArea;
	}
	public static int[] getScanListFromRT(int scanNum, double rt, double rtTolerance, String fileName) {
		MSXMLParser parser = new MSXMLParser(fileName);
		ArrayList<Integer> tmpScanList = new ArrayList<Integer>();
		double rtlb = rt - rtTolerance;
		double rtub = rt + rtTolerance;
		int prerunner = scanNum --;
		double tmprt = getRT(fileName, prerunner);
		while (tmprt > rtlb && prerunner > 0) {
			Scan scan = parser.rap(prerunner);
			ScanHeader header = scan.getHeader();
			int msLevel = header.getMsLevel();
			if (msLevel == 1) tmpScanList.add(prerunner);
			prerunner--;
			tmprt = getRT(fileName, prerunner);
		}
		int postrunner = scanNum++;
		tmprt = getRT(fileName, postrunner);
		int mxScan = returnMxScan(fileName);
		while (tmprt < rtub && postrunner < mxScan) {
			Scan scan = parser.rap(postrunner);
			ScanHeader header = scan.getHeader();
			int msLevel = header.getMsLevel();
			if (msLevel == 1) tmpScanList.add(postrunner);
			postrunner++;
			tmprt = getRT(fileName, postrunner);
		}
		int[] scanlist = new int[tmpScanList.size()];
		for (int i = 0; i < tmpScanList.size(); i++) {
			scanlist[i] = tmpScanList.get(i);
		}
		return scanlist;
	}
	public static int[] getScanList(int scanNum, String fileName) {
		MSXMLParser parser = new MSXMLParser(fileName);
		int[] scanlist = new int[40];
		//find 10 MS1 scans before the msms
		int prerunner = scanNum;
		int index = 19;
		int msLevel = -1;
		while (index >= 0) {
			while (true) {
				Scan scan = parser.rap(prerunner);
				ScanHeader header = scan.getHeader();
				msLevel = header.getMsLevel();
				if (msLevel == 1) {
					scanlist[index] = prerunner;
					prerunner--;
					break;
				}
				prerunner--;
			}
			index--;
		}
		//find 10 MS1 scans after the msms
		int postrunner = scanNum;
		index = 20;
		while (index < 40) {
			while(true) {
				Scan scan = parser.rap(postrunner);
				ScanHeader header = scan.getHeader();
				msLevel = header.getMsLevel();
				if (msLevel == 1) {
					scanlist[index] = postrunner;
					postrunner++;
					break;
				}
				postrunner++;
			}
			index++;
		}
        /*double[][] massIntensityList = scan.getMassIntensityList();
        System.out.println(msLevel);
        System.out.print("[");
        for(int i = 0; i < massIntensityList[0].length; i++)
        {
            if(i > 0)
                System.out.print(",");
            System.out.println();
            System.out.print(String.format("[%.6f,%.6f]", massIntensityList[0][i], massIntensityList[1][i]));

        }
        System.out.println("\n]");*/
        return scanlist;
	}
	public static void printScan(int scanNum, String fileName) {
		MSXMLParser parser = new MSXMLParser(fileName);
		Scan scan = parser.rap(scanNum);
		double[][] massIntensityList = scan.getMassIntensityList();
        System.out.print("[");
        for(int i = 0; i < massIntensityList[0].length; i++)
        {
            if(i > 0)
                System.out.print(",");
            System.out.println();
            System.out.print(String.format("[%.6f,%.6f]", massIntensityList[0][i], massIntensityList[1][i]));

        }
        System.out.println("\n]");
	}
	public static double integratePeakBoxCar(int[] scanlist, double mz, String fileName) {
		//int boxCarSize = 5;
		double[] peakArea = new double[scanlist.length];
		double sum = 0;
		double lB, uB;
		lB = mz - 0.000020 * mz;
		uB = mz + 0.000020 * mz;
		MSXMLParser parser = new MSXMLParser(fileName);
		int j = 0;
		for (int scanNum : scanlist) {
			Scan scan = parser.rap(scanNum);
			double[][] massIntensityList = scan.getMassIntensityList();
			for(int i = 0; i < massIntensityList[0].length; i++) {
				if (massIntensityList[0][i] < lB) continue;
				if (massIntensityList[0][i] > uB) break;
				peakArea[j] += massIntensityList[1][i];
			}
			j++;
		}
		for (int i = 0; i < peakArea.length; i++) {
			double counter = 0;
			double areaTmp = 0;
			for (int k = -2; k <= 2; k++) {
				if (i + k > 0 && i + k < peakArea.length) {
					areaTmp += peakArea[i + k];
					counter++;
				}
			}
			sum += areaTmp / counter;
		}
		return sum;
	}
}
class Precursor implements Comparable {
	public double mz;
	public double lb;
	public double ub;
	
	public Precursor() {}
	@Override
	public int compareTo(Object o) {
		if (o instanceof Precursor) {
			Precursor p = (Precursor) o;
			if (p.lb < lb) return -1;
			if (p.lb > lb) return 1;
			if (p.ub < ub) return -1;
			if (p.ub > ub) return 1;
		}
		return 0;
	}
}
