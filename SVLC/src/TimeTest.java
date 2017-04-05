import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

public class TimeTest {
	
	public static void main(String[] args) {
		String part1Path = "Henry part 1";
		String part2Path = "Henry part 2";
		
		try {
			String part1 = new String(Files.readAllBytes(Paths.get(part1Path)), StandardCharsets.UTF_8);
			String part2 = new String(Files.readAllBytes(Paths.get(part2Path)), StandardCharsets.UTF_8);
			
			part1 = part1.replace('\n', ' ').toLowerCase();
			part2 = part2.replace('\n', ' ').toLowerCase();
			
			Sequence a = Sequence.valueOf(part1);
			Sequence b = Sequence.valueOf(part2);			
			long startTime, endTime;
			
			System.out.println("Testing speed for suffixarray implementation ...");
			for (int length = 1000; length < 50000; length += 1000) {
				startTime = System.nanoTime();
				Sequence.LCSSsuffixtree(a.new Subsequence(0, length), b.new Subsequence(0, length));
				endTime = System.nanoTime();
				System.out.print(", " + (endTime - startTime));
			}
			System.out.println();
			
			System.out.println("Testing speed for suffixarray implementation ...");
			for (int length = 1000; length < 200000; length += 1000) {
				startTime = System.nanoTime();
				Sequence.LCSSdynamic(a.new Subsequence(0, length), b.new Subsequence(0, length));
				endTime = System.nanoTime();
				System.out.print(", " + (endTime - startTime));
			}
			System.out.println();
			
			
		} catch (IOException e) {
			System.err.println(e);
		}
		

		
	}

}
