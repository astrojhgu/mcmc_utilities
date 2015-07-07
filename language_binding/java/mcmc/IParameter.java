package mcmc;

public interface IParameter{
    double getValue(int i);
    void setValue(int i,double v);
    int getSize();
    void setSize(int n);
}
