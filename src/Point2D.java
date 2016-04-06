
public class Point2D{
	/*****
	*
	*   Point2D.js
	*
	*   copyright 2001-2002, Kevin Lindsey
	*
	*****/

	/*****
	*
	*   Point2D
	*
	*****/
	
	double x, y;
	
	/*****
	*
	*   constructor
	*
	*****/
	public Point2D(double x, double y) {
	        this.x = x;
	        this.y = y;
	}


	/*****
	*
	*   clone
	*
	*****/
	public Point2D clone() {
	    return new Point2D(this.x, this.y);
	};


	/*****
	*
	*   add
	*
	*****/
	public Point2D add(Point2D that) {
	    return new Point2D(this.x+that.x, this.y+that.y);
	};


	/*****
	*
	*   addEquals
	*
	*****/
	public Point2D addEquals(Point2D that) {
	    this.x += that.x;
	    this.y += that.y;

	    return this;
	};


	/*****
	*
	*   offset - used in dom_graph
	*
	*   This method is based on code written by Walter Korman
	*      http://www.go2net.com/internet/deep/1997/05/07/body.html 
	*   which is in turn based on an algorithm by Sven Moen
	*
	*****/
	public double offset(Point2D a, Point2D b) {
	    double result = 0;

	    if ( !( b.x <= this.x || this.x + a.x <= 0 ) ) {
	        double t = b.x * a.y - a.x * b.y;
	        double s;
	        double d;

	        if ( t > 0 ) {
	            if ( this.x < 0 ) {
	                s = this.x * a.y;
	                d = s / a.x - this.y;
	            } else if ( this.x > 0 ) {
	                s = this.x * b.y;
	                d = s / b.x - this.y;
	            } else {
	                d = -this.y;
	            }
	        } else {
	            if ( b.x < this.x + a.x ) {
	                s = ( b.x - this.x ) * a.y;
	                d = b.y - (this.y + s / a.x);
	            } else if ( b.x > this.x + a.x ) {
	                s = (a.x + this.x) * b.y;
	                d = s / b.x - (this.y + a.y);
	            } else {
	                d = b.y - (this.y + a.y);
	            }
	        }

	        if ( d > 0 ) {
	            result = d;
	        }
	    }

	    return result;
	};


	/*****
	*
	*   rmoveto
	*
	*****/
	public void rmoveto(double dx, double dy) {
	    this.x += dx;
	    this.y += dy;
	};


	/*****
	*
	*   scalarAdd
	*
	*****/
	public Point2D scalarAdd(double scalar) {
	    return new Point2D(this.x+scalar, this.y+scalar);
	};


	/*****
	*
	*   scalarAddEquals
	*
	*****/
	public Point2D scalarAddEquals(double scalar) {
	    this.x += scalar;
	    this.y += scalar;

	    return this;
	};


	/*****
	*
	*   subtract
	*
	*****/
	public Point2D subtract(Point2D that) {
	    return new Point2D(this.x-that.x, this.y-that.y);
	};


	/*****
	*
	*   subtractEquals
	*
	*****/
	public Point2D subtractEquals(Point2D that) {
	    this.x -= that.x;
	    this.y -= that.y;

	    return this;
	};


	/*****
	*
	*   scalarSubtract
	*
	*****/
	public Point2D scalarSubtract(double scalar) {
	    return new Point2D(this.x-scalar, this.y-scalar);
	};


	/*****
	*
	*   scalarSubtractEquals
	*
	*****/
	public Point2D scalarSubtractEquals(double scalar) {
	    this.x -= scalar;
	    this.y -= scalar;

	    return this;
	};


	/*****
	*
	*   multiply
	*
	*****/
	public Point2D multiply (double scalar) {
	    return new Point2D(this.x*scalar, this.y*scalar);
	};


	/*****
	*
	*   multiplyEquals
	*
	*****/
	public Point2D multiplyEquals (double scalar) {
	    this.x *= scalar;
	    this.y *= scalar;

	    return this;
	};


	/*****
	*
	*   divide
	*
	*****/
	public Point2D divide (double scalar) {
	    return new Point2D(this.x/scalar, this.y/scalar);
	};


	/*****
	*
	*   divideEquals
	*
	*****/
	public Point2D divideEquals (double scalar) {
	    this.x /= scalar;
	    this.y /= scalar;

	    return this;
	};


	/*****
	*
	*   comparison methods
	*
	*   these were a nice idea, but ...  It would be better to define these names
	*   in two parts so that the first part is the x comparison and the second is
	*   the y.  For example, to test p1.x < p2.x and p1.y >= p2.y, you would call
	*   p1.lt_gte(p2).  Honestly, I only did these types of comparisons in one
	*   Intersection routine, so these probably could be removed.
	*
	*****/

	/*****
	*
	*   compare
	*
	*****/
	public double compare (Point2D that) {
		double ret = this.x - that.x;
		if (ret==0)
			ret=this.y - that.y;
	  return ret;
	};


	/*****
	*
	*   eq - equal
	*
	*****/
	public boolean eq (Point2D that) {
	    return ( this.x == that.x && this.y == that.y );
	};


	/*****
	*
	*   lt - less than
	*
	*****/
	public boolean lt (Point2D that) {
	    return ( this.x < that.x && this.y < that.y );
	};


	/*****
	*
	*   lte - less than or equal
	*
	*****/
	public boolean lte (Point2D that) {
	    return ( this.x <= that.x && this.y <= that.y );
	};


	/*****
	*
	*   gt - greater than
	*
	*****/
	public boolean gt (Point2D that) {
	    return ( this.x > that.x && this.y > that.y );
	};


	/*****
	*
	*   gte - greater than or equal
	*
	*****/
	public boolean gte (Point2D that) {
	    return ( this.x >= that.x && this.y >= that.y );
	};


	/*****
	*
	*   utility methods
	*
	*****/

	/*****
	*
	*   lerp
	*
	*****/
	public Point2D lerp (Point2D that, double t) {
	    return new Point2D(
	        this.x + (that.x - this.x) * t,
	        this.y + (that.y - this.y) * t
	    );
	};


	/*****
	*
	*   distanceFrom
	*
	*****/
	public double distanceFrom (Point2D that) {
	    double dx = this.x - that.x;
	    double dy = this.y - that.y;

	    return Math.sqrt(dx*dx + dy*dy);
	};


	/*****
	*
	*   min
	*
	*****/
	public Point2D min (Point2D that) {
	    return new Point2D(
	        Math.min( this.x, that.x ),
	        Math.min( this.y, that.y )
	    );
	};


	/*****
	*
	*   max
	*
	*****/
	public Point2D max (Point2D that) {
	    return new Point2D(
	        Math.max( this.x, that.x ),
	        Math.max( this.y, that.y )
	    );
	};


	/*****
	*
	*   toString
	*
	*****/
	public String toString () {
	    return this.x + "," + this.y;
	};


	/*****
	*
	*   get/set methods
	*
	*****/

	/*****
	*
	*   setXY
	*
	*****/
	public void setXY (double x, double y) {
	    this.x = x;
	    this.y = y;
	};


	/*****
	*
	*   setFromPoint
	*
	*****/
	public void setFromPoint (Point2D that) {
	    this.x = that.x;
	    this.y = that.y;
	};


	/*****
	*
	*   swap
	*
	*****/
	public void swap (Point2D that) {
	    double x = this.x;
	    double y = this.y;

	    this.x = that.x;
	    this.y = that.y;

	    that.x = x;
	    that.y = y;
	};

}
