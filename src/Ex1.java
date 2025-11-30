/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
		/// add you code below
         /// for 2 points
        if(lx == 2){
            ans = new double[2];
            ans[1] = (yy[1] - yy[0])/(xx[1] - xx[0]);
            ans[0] =  yy[0] - xx[0]*ans[1];
            return ans;
        }
        /// /// for 3 points
        if (lx==3){
            ans=new double[3];

            double x1 = xx[0], y1 = yy[0];
            double x2 = xx[1], y2 = yy[1];
            double x3 = xx[2], y3 = yy[2];

            double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);

            double A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
            double B = (x3 * x3 * (y1 - y2) +
                    x2 * x2 * (y3 - y1) +
                    x1 * x1 * (y2 - y3)) / denom;
            double C = (x2 * x3 * (x2 - x3) * y1 +
                    x3 * x1 * (x3 - x1) * y2 +
                    x1 * x2 * (x1 - x2) * y3) / denom;

            ans[2] = A;
            ans[1] = B;
            ans[0] = C;

            return ans;
        }
		////////////////////
		}
		return ans;
	}
	/** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true iff p1 represents the same polynomial function as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;
        /// add you code below
        if (p1 == null || p2 == null) {
            return false;
        }
        int len1 = p1.length - 1;
        int len2 = p2.length - 1;
        int n = Math.max(len1 ,len2);

        for (int i=0; i<=n ;i++) {
            double x = i;
            double v1 = f(p1, x);
            double v2 = f(p2, x);

            if (Math.abs(v1 - v2) > EPS) {
                ans= false;
                break;

            }
        }

         ////////////////////
		return ans;
	}

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
	 */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
            //// add you code below
            int len = poly.length-1;
            for (int i=len; i>=0; i--){
                if (poly[i] != 0){
                    /// // handling for signs
                    if (ans.length() > 0 && poly[i] > 0) {
                        ans += " + ";
                    } else if (poly[i] < 0) {
                        ans += " ";
                    }
                    /// / handling for nums
                    if (i == 0) {
                        ans += poly[i];
                    } else if (i == 1) {
                        ans += poly[i] + "x";
                    } else {
                        ans += poly[i] + "x^" + i;
                    }
                }

            }
            ////////////////////
		}
		return ans;
	}
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        /// add you code below
        if (x1 > x2) {
            double temp = x1;
            x1 = x2;
            x2 = temp;
        }

        double g1 = f(p1, x1) - f(p2, x1);
        double xm = (x1 + x2) / 2.0;
        double gm = f(p1, xm) - f(p2, xm);

        if (Math.abs(gm) < eps) {
            return xm;
        }

        if (g1 * gm <= 0) {
            return sameValue(p1, p2, x1, xm, eps);
        } else {
            return sameValue(p1, p2, xm, x2, eps);
        }
    }
    ////////////////////

	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
        ///add you code below
        if (x1 > x2) {
            double temp = x1;
            x1 = x2;
            x2 = temp;
        }
            double h = (x2 - x1) / numberOfSegments;
            double ans = 0.0;

            double xPrev = x1;
            double yPrev = f(p, xPrev);

            for (int i = 1; i <= numberOfSegments; i++) {
                double xCurr = x1 + i * h;
                double yCurr = f(p, xCurr);

                double xd = xCurr - xPrev;
                double yd = yCurr - yPrev;

                double dis = Math.sqrt(xd * xd + yd * yd);

                ans += dis;

                xPrev = xCurr;
                yPrev = yCurr;
            }
            return ans;
        }
         ////////////////////


	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
     *
     * Implementation summary:
     * 1. If x1 > x2, the values are swapped to ensure a valid range.
     * 2. The interval is divided into 'numberOfTrapezoid' sub-intervals (minimum 1).
     * 3. For each sub-interval, compute the difference d(x) = f(p1, x) - f(p2, x).
     * 4. If d1 and d2 have the same sign, compute the regular trapezoid area.
     * 5. If a sign change occurs (the graphs cross), find the intersection point using sameValue(),
     * split the interval at that point, and compute each partial area separately.
     * 6. The absolute value of the area is accumulated.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */

    public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		double ans = 0.0;
        ///add you code below
         if (x1 > x2) {
         double temp = x1;
         x1 = x2;
         x2 = temp;
         }
         int n = Math.max(1, numberOfTrapezoid);
        double h = (x2 - x1) / n;

        for (int i =0;i<n;i++){
            double c1 = h * i + x1;
            double c2 = h * (i+1) + x1;

            double d1 = f(p1, c1) - f(p2, c1);
            double d2 = f(p1, c2) - f(p2, c2);

            // if there is a sign change, locate the crossing point and split the trapezoid
            if (d1 * d2 < 0) {
                double xr = sameValue(p1, p2, c1, c2, EPS);
                double mid = xr;
                ans += Math.abs((d1 + 0) * 0.5 * (mid - c1));
                ans += Math.abs((0 + d2) * 0.5 * (c2 - mid));
            } else {
                double trp = Math.abs((d1 + d2) * 0.5 * h);
                ans += trp;
            }

        }


         //////////////////////
		return ans;
	}


    /**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * @param p - a String representing polynomial function.
	 * @return
     *
     * The function uses two helper functions:
     * - getPowerFromTerm(String term): extracts the power (b in ax^b)
     * - getCoefFromTerm(String term): extracts the coefficient (a in ax^b)
     *
     * Steps:
     * 1) Clean whitespace, convert "-" to "+-" so we can split correctly
     * 2) Split into terms (monoms)
     * 3) Determine maximum power to know array size
     * 4) Fill coefficient values in the correct power indices
     */

    public static double[] getPolynomFromString(String p) {
        // if null then ZERO
        if (p == null) {
            return ZERO;
        }

        p = p.trim();
        if (p.equals("")) {
            return ZERO;
        }

        ///  clean spaces
        String s = p.replace(" ", "");   // "x^2 +x +1" → "x^2+x+1"

        ///  replace " - " with " +- "
        s = s.replace("-", "+-");

        /// if + at start we take down
        if (s.startsWith("+")) {
            s = s.substring(1);
        }

        /// split into monoms
        String[] parts = s.split("\\+");   // "x^2+x+1" → ["x^2","x","1"]

        /// getting max power
        int maxPower = 0;
        for (int i = 0; i < parts.length; i++) {
            String t = parts[i];
            if (!t.equals("")) {
                int power = getPowerFromTerm(t);
                if (power > maxPower) {
                    maxPower = power;
                }
            }
        }
        double[] ans = new double[maxPower + 1];

        /// filling the coefficients in their places
        for (int i = 0; i < parts.length; i++) {
            String t = parts[i];
            if (!t.equals("")) {
                int power = getPowerFromTerm(t);
                double coef = getCoefFromTerm(t);
                ans[power] += coef;
            }
        }

        return ans;
    }

    /**
     * Extracts the power (b) from a polynomial term represented as a string.
     * Example formats:
     *  "5x^3" -> returns 3
     *  "x^2"  -> returns 2
     *  "7x"   -> returns 1
     *  "4"    -> returns 0 (no x => constant term)
     *
     * @param term - a polynomial single term string (e.g. "3x^2", "-x", "5")
     * @return the power of the term as an integer
     */

    private static int getPowerFromTerm(String term) {
        if (term == null || term.equals("")) {
            return 0;
        }

        if (!term.contains("x")) {
            return 0;
        }

        if (term.contains("^")) {
            String[] splitPower = term.split("\\^");
            if (splitPower.length > 1) {
                return Integer.parseInt(splitPower[1]);
            } else {
                return 1;
            }
        }


        return 1;
    }

    /**
     * Extracts the coefficient (a) from a polynomial term represented as a string.
     * Example formats:
     *  "5x^3" -> returns 5.0
     *  "-x^2" -> returns -1.0
     *  "x"    -> returns 1.0
     *  "-x"   -> returns -1.0
     *  "7"    -> returns 7.0
     *
     * @param term - a polynomial term (e.g. "3x^2", "-x", "+4")
     * @return coefficient as double
     */

    private static double getCoefFromTerm(String term) {
        if (term == null || term.equals("")) {
            return 0.0;
        }

        if (!term.contains("x")) {

            return Double.parseDouble(term);
        }

        // מפצלים לפי x כדי לקחת את מה שלפניו
        String[] splitCoef = term.split("x", 2);
        String coefStr = splitCoef[0];

        if (coefStr.equals("") || coefStr.equals("+")) {
            return 1.0;
        }
        if (coefStr.equals("-")) {
            return -1.0;
        }

        return Double.parseDouble(coefStr);
    }

    /**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */

	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        ///add you code below
        if (p1 != null && p2 != null) {

            int n = Math.max(p1.length, p2.length);
            ans = new double[n];

            for (int i = 0; i < n; i++) {

                double a = 0;
                double b = 0;

                if (i < p1.length) {
                    a = p1[i];
                }
                if (i < p2.length) {
                    b = p2[i];
                }

                ans[i] = a + b;
            }
        }

         /////////////////// /
		return ans;
	}
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        /// add you code below
        if (p1 != null && p2 != null) {
        int n = p1.length+ p2.length-1;
        ans = new double[n];

        for(int i=0; i<p1.length;i++){
            for(int j=0;j<p2.length;j++){
                ans[i+j] = ans[i+j] + (p1[i]*p2[j]);
            }
        }
        }
         /////////////////// */
		return ans;
	}
	/**
	 * This function computes the derivative of the p0 polynomial function.
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;//
        /// add you code below/
        if(po!=null&&po.length>1){
            int len = po.length;
            ans = new double[len-1];
            for(int i=0;i< ans.length;i++){
                ans[i] = po[i+1] * (i+1);
            }
        }

         /////////////////// */
		return ans;
	}
}
