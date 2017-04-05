package general;

public class ArrayPrinter {
	
	public static <K> String print(K[][] array) {
		StringBuilder result = new StringBuilder();
		for (int i=0; i<array.length; i++) {
			String comma = "";
			for (int j=0; j<array[i].length; j++) {
				result.append(comma);
				result.append(array[i][j]);
				comma = ", ";
			}
			result.append('\n');
		}
		return result.toString();
	}
	
	public static String print(double[][] array) {
		StringBuilder result = new StringBuilder();
		for (int i=0; i<array.length; i++) {
			String comma = "";
			for (int j=0; j<array[i].length; j++) {
				result.append(comma);
				result.append(array[i][j]);
				comma = ", ";
			}
			result.append('\n');
		}
		return result.toString();
	}

}
