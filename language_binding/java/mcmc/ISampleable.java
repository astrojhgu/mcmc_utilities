package mcmc;

public interface ISampleable{
    double evalLog(IParameter p);
    double[] varRange(IParameter p,int ndim);
    double[] initPoints(IParameter p,int ndim);
}
