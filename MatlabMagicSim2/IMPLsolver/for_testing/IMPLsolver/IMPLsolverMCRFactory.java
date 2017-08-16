/*
 * MATLAB Compiler: 6.3 (R2016b)
 * Date: Mon Jul 24 03:14:49 2017
 * Arguments: "-B" "macro_default" "-W" "java:IMPLsolver,Class1" "-T" "link:lib" "-d" 
 * "/home/denise/Masterthesis/Matlab/MatlabMagicSim2/IMPLsolver/for_testing" 
 * "class{Class1:/home/denise/Masterthesis/Matlab/MatlabMagicSim2/IMPLsolver.m}" 
 */

package IMPLsolver;

import com.mathworks.toolbox.javabuilder.*;
import com.mathworks.toolbox.javabuilder.internal.*;

/**
 * <i>INTERNAL USE ONLY</i>
 */
public class IMPLsolverMCRFactory
{
   
    
    /** Component's uuid */
    private static final String sComponentId = "IMPLsolver_F9DFBD1B57319A93B112BC16122528C4";
    
    /** Component name */
    private static final String sComponentName = "IMPLsolver";
    
   
    /** Pointer to default component options */
    private static final MWComponentOptions sDefaultComponentOptions = 
        new MWComponentOptions(
            MWCtfExtractLocation.EXTRACT_TO_CACHE, 
            new MWCtfClassLoaderSource(IMPLsolverMCRFactory.class)
        );
    
    
    private IMPLsolverMCRFactory()
    {
        // Never called.
    }
    
    public static MWMCR newInstance(MWComponentOptions componentOptions) throws MWException
    {
        if (null == componentOptions.getCtfSource()) {
            componentOptions = new MWComponentOptions(componentOptions);
            componentOptions.setCtfSource(sDefaultComponentOptions.getCtfSource());
        }
        return MWMCR.newInstance(
            componentOptions, 
            IMPLsolverMCRFactory.class, 
            sComponentName, 
            sComponentId,
            new int[]{9,1,0}
        );
    }
    
    public static MWMCR newInstance() throws MWException
    {
        return newInstance(sDefaultComponentOptions);
    }
}
