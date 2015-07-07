import mcmc.*;

class Param
    implements IParameter
{
    private double[] p=new double[2];
    public double getValue(int i)
    {
	return p[i];
    }

    public void setValue(int i,double v)
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
		result+=p.getValue(i)*p.getValue(i);
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
	Param p=new Param();
	
	for(int i=0;i<10000;++i)
	    {
		Sampler.gibbsSample(new Gdist(),p,true);
		System.out.println(""+p.getValue(0)+" "+p.getValue(1));
	    }
    }    
}

    
