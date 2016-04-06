import java.util.ArrayList;

public class Polynomial{
	/*****
	*
	*   Polynomial.js
	*
	*   copyright 2002, Kevin Lindsey
	*
	*****/
	
	private static final double LN10 = 2.30258509299;
	private static final double LN2 = 0.69314718056;
			
	static final double TOLERANCE = 1e-6;
	static final double ACCURACY  = 6;
	
	//TODO: can this be final?
	double[] coefs;

  String _variable;
  double _s;


	/*****
	*
	*   interpolate - class method
	*
	*****/
	static Derivative interpolate(double[] xs, double[] ys, int n, int offset, double x) {
	    //if ( xs.constructor !== Array || ys.constructor !== Array )
	        //throw new Error("Polynomial.interpolate: xs and ys must be arrays");
	    if ( Double.isNaN(n) || Double.isNaN(offset) || Double.isNaN(x) )
	        throw new Error("Polynomial.interpolate: n, offset, and x must be numbers");

	    double y  = 0;
	    double dy = 0;
	    double[] c = new double[n];
	    double[] d = new double[n];
	    int ns = 0;
	    Derivative result;

	    double diff = Math.abs(x - xs[offset]);
	    for ( int i = 0; i < n; i++ ) {
	        double dift = Math.abs(x - xs[offset+i]);

	        if ( dift < diff ) {
	            ns = i;
	            diff = dift;
	        }
	        c[i] = d[i] = ys[offset+i];
	    }
	    y = ys[offset+ns];
	    ns--;

	    for ( int m = 1; m < n; m++ ) {
	        for ( int i = 0; i < n-m; i++ ) {
	        	double ho = xs[offset+i] - x;
	        	double hp = xs[offset+i+m] - x;
	        	double w = c[i+1]-d[i];
	        	double den = ho - hp;

	            if ( den == 0.0 ) {
	                result = new Derivative(0, 0);
	                break;
	            }

	            den = w / den;
	            d[i] = hp*den;
	            c[i] = ho*den;
	        }
	        dy = (2*(ns+1) < (n-m)) ? c[ns+1] : d[ns--];
	        y += dy;
	    }

	    return new Derivative(y, dy);
	}


	/*****
	*
	*   constructor
	*
	*****/
	Polynomial(double... arguments) {
	    this.init( arguments );
	}


	/*****
	*
	*   init
	*
	*****/
	public void init(double... coefs) {
	    this.coefs = new double[coefs.length];

	    for ( int i = coefs.length - 1; i >= 0; i-- )
	        this.coefs[coefs.length - 1 - i]=coefs[i];

	    this._variable = "t";
	    this._s = 0;
	}


	/*****
	*
	*   eval
	*
	*****/
	public double eval(double x) {
	    if ( Double.isNaN(x) )
	        throw new Error("Polynomial.eval: parameter must be a number");

	    double result = 0;

	    for ( int i = this.coefs.length - 1; i >= 0; i-- )
	        result = result * x + this.coefs[i];

	    return result;
	}


	/*****
	*
	*   add
	*
	*****/
	public Polynomial add(Polynomial that) {
	    Polynomial result = new Polynomial();
	    int d1 = this.getDegree();
	    int d2 = that.getDegree();
	    int dmax = Math.max(d1,d2);

	    for ( int i = 0; i <= dmax; i++ ) {
	        double v1 = (i <= d1) ? this.coefs[i] : 0;
	        double v2 = (i <= d2) ? that.coefs[i] : 0;

	        result.coefs[i] = v1 + v2;
	    }

	    return result;
	}


	/*****
	*
	*   multiply
	*
	*****/
	public Polynomial multiply(Polynomial that) {
		Polynomial result = new Polynomial();

	        result.coefs = new double[this.getDegree() + that.getDegree()];

	    for ( int i = 0; i <= this.getDegree(); i++ )
	        for ( int j = 0; j <= that.getDegree(); j++ )
	            result.coefs[i+j] += this.coefs[i] * that.coefs[j];

	    return result;
	}


	/*****
	*
	*   divide_scalar
	*
	*****/
	public void divide_scalar(double scalar) {
	    for ( int i = 0; i < this.coefs.length; i++ )
	        this.coefs[i] /= scalar;
	}


	/*****
	*
	*   simplify
	*
	*****/
//	Polynomial.prototype.simplify = function() {
//    for ( var i = this.getDegree(); i >= 0; i-- ) {
//        if ( Math.abs( this.coefs[i] ) <= Polynomial.TOLERANCE )
//            this.coefs.pop();
//        else
//            break;
//    }
//};
	public void simplify() {
	    for ( int i = this.getDegree(); i >= 0; i-- ) {
	        if ( Math.abs( this.coefs[i] ) <= Polynomial.TOLERANCE ){
	        	double[] old = this.coefs;
	        	coefs = new double[coefs.length-1];
	        	for (int j=0; j<coefs.length; j++)
	        		coefs[j] = old[j];
	        }
	        else
	            break;
	    }
	};


	/*****
	*
	*   bisection
	*
	*****/
	public Double bisection(double min, double max) {
	    double minValue = this.eval(min);
	    double maxValue = this.eval(max);
	    Double result = null;
	    
	    if ( Math.abs(minValue) <= Polynomial.TOLERANCE )
	        result = min;
	    else if ( Math.abs(maxValue) <= Polynomial.TOLERANCE )
	        result = max;
	    else if ( minValue * maxValue <= 0 ) {
	        double tmp1  = Math.log(max - min);
	        double tmp2  = LN10 * Polynomial.ACCURACY;
	        double iters = Math.ceil( (tmp1+tmp2) / LN2 );

	        for ( int i = 0; i < iters; i++ ) {
	            result = 0.5 * (min + max);
	            double value = this.eval(result);

	            if ( Math.abs(value) <= Polynomial.TOLERANCE ) {
	                break;
	            }

	            if ( value * minValue < 0 ) {
	                max = result;
	                maxValue = value;
	            } else {
	                min = result;
	                minValue = value;
	            }
	        }
	    }
	    return result;
	};


	/*****
	*
	*   toString
	*
	*****/
	public String toString() {
	    String[] coefs = new String[this.coefs.length];
	    String[] signs = new String[this.coefs.length];
	    
	    for ( int i = this.coefs.length - 1; i >= 0; i-- ) {
	        double value = Math.round(this.coefs[i]*1000)/1000.0;
	        //var value = this.coefs[i];

	        if ( value != 0 ) {
	            String sign = ( value < 0 ) ? " - " : " + ";

	            value = Math.abs(value);
	            String value2 = "";
	            if ( i > 0 )
	                if ( value == 1 )
	                    value2 = this._variable;
	                else
	                    value2 += this._variable;
	            if ( i > 1 )
	            	value2 += "^" + i;

	            signs[i]=sign;
	            coefs[i]=value2;
	        }
	    }

	    signs[0] = signs[0].equals(" + ") ? "" : "-";

	    String result = "";
	    for ( int i = 0; i < coefs.length; i++ )
	        result += signs[i] + coefs[i];
	    
	    return result;
	};


	/*****
	*
	*   trapezoid
	*   Based on trapzd in "Numerical Recipes in C", page 137
	*
	*****/
	public double trapezoid(double min, double max, int n) {
	    if ( Double.isNaN(min) || Double.isNaN(max) || Double.isNaN(n) )
	        throw new Error("Polynomial.trapezoid: parameters must be numbers");

	    double range = max - min;
	    double TOLERANCE = 1e-7;

	    if ( n == 1 ) {
	    	double minValue = this.eval(min);
	    	double maxValue = this.eval(max);
	        this._s = 0.5*range*( minValue + maxValue );
	    } else {
	        int it = 1 << (n-2);
	        double delta = range / it;
	        double x = min + 0.5*delta;
	        double sum = 0;

	        for ( int i = 0; i < it; i++ ) {
	            sum += this.eval(x);
	            x += delta;
	        }
	        this._s = 0.5*(this._s + range*sum/it);
	    }

	    if ( Double.isNaN(this._s) )
	        throw new Error("Polynomial.trapezoid: this._s is NaN");

	    return this._s;
	};


	/*****
	*
	*   simpson
	*   Based on trapzd in "Numerical Recipes in C", page 139
	*
	*****/
	public double simpson(double min, double max) {
	    if ( Double.isNaN(min) || Double.isNaN(max) )
	        throw new Error("Polynomial.simpson: parameters must be numbers");

	    double range = max - min;
	    double st = 0.5 * range * ( this.eval(min) + this.eval(max) );
	    double t = st;
	    double s = 4.0*st/3.0;
	    double os = s;
	    double ost = st;
	    double TOLERANCE = 1e-7;

	    int it = 1;
	    for ( int n = 2; n <= 20; n++ ) {
	        double delta = range / it;
	        double x     = min + 0.5*delta;
	        double sum   = 0;

	        for ( int i = 1; i <= it; i++ ) {
	            sum += this.eval(x);
	            x += delta;
	        }

	        t = 0.5 * (t + range * sum / it);
	        st = t;
	        s = (4.0*st - ost)/3.0;

	        if ( Math.abs(s-os) < TOLERANCE*Math.abs(os) )
	            break;

	        os = s;
	        ost = st;
	        it <<= 1;
	    }

	    return s;
	};


	/*****
	*
	*   romberg
	*
	*****/
	public double romberg(double min, double max) {
	    if ( Double.isNaN(min) || Double.isNaN(max) )
	        throw new Error("Polynomial.romberg: parameters must be numbers");

	    int MAX = 20;
	    int K = 3;
	    double TOLERANCE = 1e-6;
	    double[] s = new double[MAX+1];
	    double[] h = new double[MAX+1];
	    Derivative result = new Derivative(0, 0);

	    h[0] = 1.0;
	    for ( int j = 1; j <= MAX; j++ ) {
	        s[j-1] = this.trapezoid(min, max, j);
	        if ( j >= K ) {
	            result = Polynomial.interpolate(h, s, K, j-K, 0.0);
	            if ( Math.abs(result.dy) <= TOLERANCE*result.y) break;
	        }
	        s[j] = s[j-1];
	        h[j] = 0.25 * h[j-1];
	    }

	    return result.y;
	};


	/*****
	*
	*   get/set methods
	*
	*****/

	/*****
	*
	*   get degree
	*
	*****/
	public int getDegree() {
	    return this.coefs.length - 1;
	};


	/*****
	*
	*   getDerivative
	*
	*****/
	public Polynomial getDerivative() {
	    Polynomial derivative = new Polynomial();

	    derivative.coefs = new double[this.coefs.length-1];
	    for ( int i = 1; i < this.coefs.length; i++ ) {
	        derivative.coefs[i-1] = (i*this.coefs[i]);
	    }

	    return derivative;
	};


	/*****
	*
	*   getRoots
	*
	*****/
	public double[] getRoots() {
	    double[] result;

	    this.simplify();
	    switch ( this.getDegree() ) {
	        case 0: result = new double[0];            break;
	        case 1: result = this.getLinearRoot();     break;
	        case 2: result = this.getQuadraticRoots(); break;
	        case 3: result = this.getCubicRoots();     break;
	        case 4: result = this.getQuarticRoots();   break;
	        default:
	            result = new double[0];
	            // should try Newton's method and/or bisection
	    }

	    return result;
	};


	/*****
	*
	*   getRootsInInterval
	*
	*****/
	public double[] getRootsInInterval(double min, double max) {
	    ArrayList<Double> roots = new ArrayList<Double>();

	    if ( this.getDegree() == 1 ) {
	        Double root = this.bisection(min, max);
	        if ( root != null )
	        	roots.add(root);
	    } else {
	        // get roots of derivative
	        Polynomial deriv  = this.getDerivative();
	        double[] droots = deriv.getRootsInInterval(min, max);
	        

	        if ( droots.length > 0 ) {
	            // find root on [min, droots[0]]
	            Double root = this.bisection(min, droots[0]);
	            if ( root != null )
	            	roots.add(root);

	            // find root on [droots[i],droots[i+1]] for 0 <= i <= count-2
	            for ( int i = 0; i <= droots.length-2; i++ ) {
	                root = this.bisection(droots[i], droots[i+1]);
	                if ( root != null )
	                	roots.add(root);
	            }

	            // find root on [droots[count-1],xmax]
	            root = this.bisection(droots[droots.length-1], max);
	            if ( root != null )
	            	roots.add(root);
	        } else {
	            // polynomial is monotone on [min,max], has at most one root
	            Double root = this.bisection(min, max);
	            if ( root != null )
	            	roots.add(root);
	        }
	    }
	    
	    double[] ret = new double[roots.size()];
	    for (int i=0; i<roots.size(); i++)
	    	ret[i] = roots.get(i);
	    return ret;
	};


	/*****
	*
	*   getLinearRoot
	*
	*****/
	public double[] getLinearRoot() {
	    double[] result;
	    double a = this.coefs[1];
	    
	    if ( a != 0 )
	    	result = new double[]{( -this.coefs[0] / a )};
	    else
	    	result = new double[0];

	    return result;
	};


	/*****
	*
	*   getQuadraticRoots
	*
	*****/
	public double[] getQuadraticRoots() {
	    double[] results = new double[0];

	    if ( this.getDegree() == 2 ) {
	        double a = this.coefs[2];
	        double b = this.coefs[1] / a;
	        double c = this.coefs[0] / a;
	        double d = b*b - 4*c;

	        if ( d > 0 ) {
	        	double e = Math.sqrt(d);
	            
	            results = new double[]{
	            		( 0.5 * (-b + e) ),
	            		( 0.5 * (-b - e) )};
	        } else if ( d == 0 ) {
	            // really two roots with same value, but we only return one
	            results = new double[]{( 0.5 * -b )};
	        }
	    }

	    return results;
	};


	/*****
	*
	*   getCubicRoots
	*
	*   This code is based on MgcPolynomial.cpp written by David Eberly.  His
	*   code along with many other excellent examples are avaiable at his site:
	*   http://www.magic-software.com
	*
	*****/
	public double[] getCubicRoots() {
	    double[] results = new double[0];

	    if ( this.getDegree() == 3 ) {
	        double c3 = this.coefs[3];
	        double c2 = this.coefs[2] / c3;
	        double c1 = this.coefs[1] / c3;
	        double c0 = this.coefs[0] / c3;

	        double a       = (3*c1 - c2*c2) / 3;
	        double b       = (2*c2*c2*c2 - 9*c1*c2 + 27*c0) / 27;
	        double offset  = c2 / 3;
	        double discrim = b*b/4 + a*a*a/27;
	        double halfB   = b / 2;

	        if ( Math.abs(discrim) <= Polynomial.TOLERANCE ) 
	        	discrim = 0;
	        
	        if ( discrim > 0 ) {
	        	double e = Math.sqrt(discrim);
	            double tmp;
	            double root;

	            tmp = -halfB + e;
	            if ( tmp >= 0 )
	                root = Math.pow(tmp, 1/3.0);
	            else
	                root = -Math.pow(-tmp, 1/3.0);

	            tmp = -halfB - e;
	            if ( tmp >= 0 )
	                root += Math.pow(tmp, 1/3.0);
	            else
	                root -= Math.pow(-tmp, 1/3.0);

	            results = new double[]{ root - offset };
	        } else if ( discrim < 0 ) {
	            double distance = Math.sqrt(-a/3);
	            double angle    = Math.atan2( Math.sqrt(-discrim), -halfB) / 3;
	            double cos      = Math.cos(angle);
	            double sin      = Math.sin(angle);
	            double sqrt3    = Math.sqrt(3);

	            results = new double[]{
	            		( 2*distance*cos - offset ),
	            		( -distance * (cos + sqrt3 * sin) - offset),
	            		( -distance * (cos - sqrt3 * sin) - offset)};
	        } else {
	            double tmp;

	            if ( halfB >= 0 )
	                tmp = -Math.pow(halfB, 1/3.0);
	            else
	                tmp = Math.pow(-halfB, 1/3.0);

	            results = new double[]{
	            		( 2*tmp - offset ),
	            		// really should return next root twice, but we return only one
	            		( -tmp - offset )};
	        }
	    }

	    return results;
	};


	/*****
	*
	*   getQuarticRoots
	*
	*   This code is based on MgcPolynomial.cpp written by David Eberly.  His
	*   code along with many other excellent examples are avaiable at his site:
	*   http://www.magic-software.com
	*
	*****/
	public double[] getQuarticRoots() {
	    double[] results = new double[0];

	    if ( this.getDegree() == 4 ) {
	        double c4 = this.coefs[4];
	        double c3 = this.coefs[3] / c4;
	        double c2 = this.coefs[2] / c4;
	        double c1 = this.coefs[1] / c4;
	        double c0 = this.coefs[0] / c4;

	        double[] resolveRoots = new Polynomial(
	            1, -c2, c3*c1 - 4*c0, -c3*c3*c0 + 4*c2*c0 -c1*c1
	        ).getCubicRoots();
	        double y       = resolveRoots[0];
	        double discrim = c3*c3/4 - c2 + y;

	        if ( Math.abs(discrim) <= Polynomial.TOLERANCE ) discrim = 0;

	        if ( discrim > 0 ) {
	        	double e     = Math.sqrt(discrim);
	        	double t1    = 3*c3*c3/4 - e*e - 2*c2;
	        	double t2    = ( 4*c3*c2 - 8*c1 - c3*c3*c3 ) / ( 4*e );
	        	double plus  = t1+t2;
	        	double minus = t1-t2;

	            if ( Math.abs(plus)  <= Polynomial.TOLERANCE ) plus  = 0;
	            if ( Math.abs(minus) <= Polynomial.TOLERANCE ) minus = 0;

	            if ( plus >= 0 ) {
	            	double f = Math.sqrt(plus);

	                results = new double[]{
	                	( -c3/4 + (e+f)/2 ),
	                	( -c3/4 + (e-f)/2 )};
	            }
	            if ( minus >= 0 ) {
	                double f = Math.sqrt(minus);

	                results = new double[]{
	                		( -c3/4 + (f-e)/2 ),
	                		( -c3/4 - (f+e)/2 )};
	            }
	        } else if ( discrim < 0 ) {
	            // no roots
	        } else {
	            double t2 = y*y - 4*c0;

	            if ( t2 >= -Polynomial.TOLERANCE ) {
	                if ( t2 < 0 ) t2 = 0;

	                t2 = 2*Math.sqrt(t2);
	                double t1 = 3*c3*c3/4 - 2*c2;
	                if ( t1+t2 >= Polynomial.TOLERANCE ) {
	                    double d = Math.sqrt(t1+t2);

	                    results = new double[]{
	                    		( -c3/4 + d/2 ),
	                    		( -c3/4 - d/2 )};
	                }
	                if ( t1-t2 >= Polynomial.TOLERANCE ) {
	                    double d = Math.sqrt(t1-t2);

	                    results = new double[]{
	                    		( -c3/4 + d/2 ),
	                    		( -c3/4 - d/2 )};
	                }
	            }
	        }
	    }

	    return results;
	};
}
