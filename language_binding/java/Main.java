import mcmc.*;

class Param
    implements IParameter
{
    private double[] p=new double[2];
    public double get(int i)
    {
	return p[i];
    }

    public void set(int i,double v)
    {
	p[i]=v;
    }

    public int getSize()
    {
	return 2;
    }

    public void setSize(int n)
    {
	p=new double[n];
    }
}

class Gdist
    implements ISampleable
{
    public double evalLog(IParameter p)
    {
	double result=0;
	for(int i=0;i<2;++i)
	    {
		result+=p.get(i)*p.get(i);
	    }
	return -result;
    }

    public double[] varRange(IParameter p,int ndim)
    {
	double[] result=new double[2];
	result[0]=-10;
	result[1]=10;
	return result;
    }

    public double[] initPoints(IParameter p,int ndim)
    {
	double[] result={-5,-1,0,1,5};
	return result;
    }
}


public class Main{
    public static void main(String[] argv)
    {
	System.setProperty("java.library.path",System.getProperty("java.library.path")+":/home/astrojhgu/src/mcmc_utilities/language_binding/java/");
	Param p=new Param();
	
	for(int i=0;i<10000;++i)
	    {
		Sampler.gibbsSample(new Gdist(),p,true);
		System.out.println(""+p.get(0)+" "+p.get(1));
	    }
    }    
}

    
