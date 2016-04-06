
public class IntersectionTest{
	
	public static void main(String[] args){

		//				Point2D a1 = new Point2D(0, 0);
		//				Point2D a2 = new Point2D(0.4, 0.4);
		//				Point2D a3 = new Point2D(0.6, 0.6);
		//				Point2D a4 = new Point2D(1, 1);
		//		
		//				Point2D b1 = new Point2D(0, 1);
		//				Point2D b2 = new Point2D(0.4, 0.6);
		//				Point2D b3 = new Point2D(0.6, 0.4);
		//				Point2D b4 = new Point2D(1, 0);

		Point2D a1 = new Point2D(0, 0);
		Point2D a2 = new Point2D(1, 0);
		Point2D a3 = new Point2D(1, 1);
		Point2D a4 = new Point2D(0, 1);

		Point2D b1 = new Point2D(1, 0);
		Point2D b2 = new Point2D(1, 1);
		Point2D b3 = new Point2D(0, 1);
		Point2D b4 = new Point2D(0, 0);

		//		Point2D a1 = new Point2D(0, 0);
		//		Point2D a2 = new Point2D(1, 0);
		//		Point2D a3 = new Point2D(2, 0);
		//		Point2D a4 = new Point2D(3, 0);
		//
		//		Point2D b1 = new Point2D(0, 0);
		//		Point2D b2 = new Point2D(0, 1);
		//		Point2D b3 = new Point2D(0, 2);
		//		Point2D b4 = new Point2D(0, 3);

		long time = System.currentTimeMillis();

		Intersection inter = Intersection.intersectBezier3Bezier3(a1, a2, a3, a4, b1, b2, b3, b4);
		for (int i = 0; i < 10000; i++)
			inter = Intersection.intersectBezier3Bezier3(a1, a2, a3, a4, b1, b2, b3, b4);

		time = System.currentTimeMillis() - time;

		System.out.println("time: " + time + "ms, " + (time / 10) + "microseconds per intersection");

		System.out.println(inter.points + " " + inter.status);
	}
	
	public void testBezier3Bezier3(){
		
	}
	
}
