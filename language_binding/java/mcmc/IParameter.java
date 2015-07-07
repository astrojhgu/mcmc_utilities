package mcmc;

public interface IParameter{
    double get(int i);
    void set(int i,double v);
    int getSize();
    void setSize(int n);
}
