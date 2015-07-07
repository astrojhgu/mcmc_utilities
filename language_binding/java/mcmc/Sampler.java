package mcmc;

public class Sampler{
    static {
	try {
	    System.loadLibrary("mcmc");
	} catch (UnsatisfiedLinkError e) {
	    System.err.println("Native code library failed to load. \n" + e);
	    System.err.println(System.getProperty("java.library.path"));
	    System.exit(1);
	}
    }
    
    public native static void gibbsSample(ISampleable d,IParameter p,boolean doMetrop);
}
