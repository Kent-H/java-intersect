
public class Vector2D{
	/*****
	*
	*   Vector2D.js
	*
	*   copyright 2001-2002, Kevin Lindsey
	*
	*****/
	
	double x, y;

	/*****
	*
	*   constructor
	*
	*****/
	Vector2D(double x, double y) {
	        this.x = x;
	        this.y = y;
	}


	/*****
	*
	*   length
	*
	*****/
	public double length() {
	    return Math.sqrt(this.x*this.x + this.y*this.y);
	};


	/*****
	*
	*   dot
	*
	*****/
	public double dot(Vector2D that) {
	    return this.x*that.x + this.y*that.y;
	};


	/*****
	*
	*   cross
	*
	*****/
	public double cross(Vector2D that) {
	    return this.x*that.y - this.y*that.x;
	}


	/*****
	*
	*   unit
	*
	*****/
	public Vector2D unit() {
	    return this.divide( this.length() );
	};


	/*****
	*
	*   unitEquals
	*
	*****/
	public Vector2D unitEquals() {
	    this.divideEquals( this.length() );

	    return this;
	};


	/*****
	*
	*   add
	*
	*****/
	public Vector2D add(Vector2D that) {
	    return new Vector2D(this.x + that.x, this.y + that.y);
	};


	/*****
	*
	*   addEquals
	*
	*****/
	public Vector2D addEquals(Vector2D that) {
	    this.x += that.x;
	    this.y += that.y;

	    return this;
	};


	/*****
	*
	*   subtract
	*
	*****/
	public Vector2D subtract(Vector2D that) {
	    return new Vector2D(this.x - that.x, this.y - that.y);
	};


	/*****
	*
	*   subtractEquals
	*
	*****/
	public Vector2D subtractEquals(Vector2D that) {
	    this.x -= that.x;
	    this.y -= that.y;

	    return this;
	};


	/*****
	*
	*   multiply
	*
	*****/
	public Vector2D multiply(double scalar) {
	    return new Vector2D(this.x * scalar, this.y * scalar);
	};


	/*****
	*
	*   multiplyEquals
	*
	*****/
	public Vector2D multiplyEquals(double scalar) {
	    this.x *= scalar;
	    this.y *= scalar;

	    return this;
	};


	/*****
	*
	*   divide
	*
	*****/
	public Vector2D divide(double scalar) {
	    return new Vector2D(this.x / scalar, this.y / scalar);
	};


	/*****
	*
	*   divideEquals
	*
	*****/
	public Vector2D divideEquals(double scalar) {
	    this.x /= scalar;
	    this.y /= scalar;

	    return this;
	};


	/*****
	*
	*   perp
	*
	*****/
	public Vector2D perp() {
	    return new Vector2D(-this.y, this.x);
	};


	/*****
	*
	*   perpendicular
	*
	*****/
	public Vector2D perpendicular(Vector2D that) {
	    return this.subtract(this.project(that));
	};


	/*****
	*
	*   project
	*
	*****/
	public Vector2D project(Vector2D that) {
	    double percent = this.dot(that) / that.dot(that);

	    return that.multiply(percent);
	};


	/*****
	*
	*   toString
	*
	*****/
	public String toString() {
	    return this.x + "," + this.y;
	};


	/*****
	*
	*   fromPoints
	*
	*****/
	public static Vector2D fromPoints(Point2D p1, Point2D p2) {
	    return new Vector2D(
	        p2.x - p1.x,
	        p2.y - p1.y
	    );
	};
}
