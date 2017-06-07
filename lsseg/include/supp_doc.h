//===========================================================================
// The Level-Set Segmentation Library (LSSEG)
//
//
// Copyright (C) 2000-2005 SINTEF ICT, Applied Mathematics, Norway.
//
// This program is free software; you can redistribute it and/or          
// modify it under the terms of the GNU General Public License            
// as published by the Free Software Foundation version 2 of the License. 
//
// This program is distributed in the hope that it will be useful,        
// but WITHOUT ANY WARRANTY; without even the implied warranty of         
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
// GNU General Public License for more details.                           
//
// You should have received a copy of the GNU General Public License      
// along with this program; if not, write to the Free Software            
// Foundation, Inc.,                                                      
// 59 Temple Place - Suite 330,                                           
// Boston, MA  02111-1307, USA.                                           
//
// Contact information: e-mail: tor.dokken@sintef.no                      
// SINTEF ICT, Department of Applied Mathematics,                         
// P.O. Box 124 Blindern,                                                 
// 0314 Oslo, Norway.                                                     
// 
//
// Other licenses are also available for this software, notably licenses
// for:
// - Building commercial software.                                        
// - Building software whose source code you wish to keep private.        
//
//===========================================================================
//===========================================================================
//                                                                           
// File: supp_doc.h                                                          
//                                                                           
// Created: Wed Sep  6 14:11:46 2006                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: supp_doc.h,v 1.27 2006/11/25 20:23:19 oan Exp $
//                                                                           
/// \file
/// \brief File containing supplementary documentation for Doxygen.
///        This file is a header file only to trigger the Doxygen 
///        documentation generator.  It is NOT to be included anywhere.
//                                                                           
//===========================================================================

//============================MAINPAGE=======================================
/// \mainpage The Level-Set Segmentation Library (LSSEG) 
/// <CENTER> <H2> With Anisotropic Smoothing Extension </H2> </CENTER>
/// <CENTER> <I>Written by Odd Andersen </I> </CENTER>
/// \section section_intro What is it?  What does it do?
/// 
/// The LSSEG library was written during 2005/2006 in the context of a shared
/// <a href="http://www.sintef.no/content/page1____8999.aspx">internal study project</a>
/// on PDE-base image processing, between the two 
/// <a href="http://www.sintef.no/content/page2____334.aspx">SINTEF ICT</a> departments 
/// <a href="http://www.sintef.no/content/page3____342.aspx">Applied Mathematics</a> 
/// and <a href="http://www.sintef.no/content/page3____349.aspx">Optical Measurement 
/// Systems and Data Analysis</a>.
///
/// The main two applications considered was <em> anisotropic image smoothing </em>
/// and <em>automatic, level-set based image segmentation </em>.  The LSSEG library
/// contains algorithms and data structures that are useful when working in these
/// contexts, and it also contains a \ref sample_programs_page "set of sample programs"
///  where it demonstrates what can be done using the tools it provides.
///
/// \section section_howtoread How to read this documentation
///
/// This documentation is put together for the sole purpose of making it easier for 
/// the unfamiliar programmer to start using the contents of this library.  It is \em not
/// meant to be used as a theoretical introduction to PDE-based image processing.
/// We have, however, provided links to \ref reference_page "a few hand-picked publications" 
/// that explain the underlying theory of what we have done.  We benefitted a lot from these 
/// publications when putting together this library, and want to express our thanks to 
/// the respective authors.
///
/// In order to get the general "feel" of this software's capabilities, we would suggest:
///
/// \li read the pages referred to by 
/// \ref section_segmentation_introductory "How Segmentation in LSSEG Works" and
/// \ref section_smoothing_introductory "How Anisotropic Smoothing in LSSEG Works".
/// \li compile and run some of the provided \ref sample_programs_page "sample programs",
///     then have a look at their code and see how they make use of the library.
/// \li study the documentation of the \ref lsseg "main classes and algorithms", in 
///     particular the \ref lsseg::Image "Image" class, and (for segmentation) the
///     \ref lsseg::LevelSetFunction "LevelSetFunction" and \ref lsseg::ForceGenerator "ForceGenerator" 
///     classes.  Some relevant algorithms would be \ref anisotropic_smoothing(), 
///     \ref compute_structure_tensor_2D() "compute_structure_tensor_XD()",
///     \ref compute_smoothing_geometry_2D() "compute_smoothing_geometry_XD()",
///     \ref LIC_2D_FS() "LIC_XD_FS()" (related to smoothing) and 
///     \ref normal_direction_flow_2D() "normal_direction_flow_XD()", 
///     \ref develop_single_region_2D() "develop_single_region_XD()" and
///     \ref develop_multiregion_2D() "develop_multiregion_2D()" (related to level-set functions and
///     image segmentation).
/// \li Have a good look at the referenced publications, particularly those chapters and sections
///     mentioned elsewhere in this documentation.  Nothing helps more than having a good theoretical 
///     understanding about what is going on.
///
/// \note This documentation is \em not written to be of a quality
/// comparable to a scholarly paper.  The list of references is \em not extensive.
/// Rather than aiming to present the works in which various ideas were 
/// \em first developed, its purpose is to provide the reader with a limited
/// reading list that would allow him/her to get a good grasp of the main
/// underlying theoretical ideas upon which this software has been developed.
///
/// \section section_segmentation_introductory How Segmentation in LSSEG Works
///
/// We have written \ref page_HowSegmentationWorks "an introduction on the principles of segmentation in LSSEG", 
/// explaining the main concepts and algorithms.  You can find it \ref page_HowSegmentationWorks "here".
///
/// \section section_smoothing_introductory How Anisotropic Smoothing in LSSEG Works
///
/// We have also written \ref page_HowSmoothingWorks "an introduction to the algorithm used for anisotropic smoothing"
/// in LSSEG.  It is to be read \ref page_HowSmoothingWorks "here".
///
/// \section section_sample_programs Sample programs
/// The sample programs found on \ref sample_programs_page "this page" make use of the library, and 
/// may be a natural starting point for exploring its capabilities and structure.  
///
/// \section additional_help Other helpful parts of this documentation
///
/// \li We have assembled a small \ref limited_glossary "glossary" of terms used 
///     in this documentation.
///
/// \li A page of \ref reference_page "references" to works consulted can be found here.
/// 
/// 
//====================NAMESPACE DOCUMENTATIONS===============================
/// \namespace std
/// \brief Namespace for the C++ Standard Library.

/// \namespace boost 
/// \brief Namespace of the Boost C++ Libraries (<I> external - not part of our project</I>)
///
/// Boost provides free peer-reviewed portable C++ source libraries, downloadable
/// from <a href="http://boost.org">http://boost.org</a>.

/// \namespace boost::lambda 
/// \brief Namespace of the Boost.Lambda library (<I> external - not part of our project</I>)
/// 
/// The Boost.Lambda library is contained within the \ref boost "Boost" C++ Libraries.
/// It implements form of <I> lambda abstractions </I> for C++.  The term originates
/// from functional programming and lambda calculus, where a lambda abstraction defines
/// an unnamed function.  The primary motivation for the BLL is to provide flexible
/// and convenient means to define unnamed function objects for STL algorithms.  Read
/// more at <a href="http://www.boost.org/doc/html/lambda.html">
/// http://www.boost.org/doc/html/lambda.html</a>.

/// \namespace cimg_library
/// \brief Namespace of the CImg library (<I> external - not part of our project</I>)
///
/// The CImg library is a C++ library for image manipulation.  It was written by 
/// David Tschumperle, and can be downloaded <a href="http://cimg.sourceforge.net/">here</a>.
/// Our project \ref lsseg mainly depends on this library for image loading and visualization
/// purposes.  
/// \note All of \em our project code that depends on CImg in one way or another is grouped
/// together in the files \ref cimg_dependent.h (declarations) and \ref cimg_dependent.C (definintions). 
/// Independence from the CImg library can be achieved by re-implementing these functions.

/// \namespace lsseg
/// \brief This is <b> the main namespace </b> of our project.  
///
/// The name \c lsseg stands for "level-set segmentation", which was the 
/// initial aim of this code project.  However, in addition to level-set 
/// segmentation code, based on \ref anchor_Brox05 "[Brox05]", this namespace also contains
/// functionality for anisotropic smoothing of images, based on an algorithm
/// presented by David Tschumperl&eacute; in \ref anchor_Tschumperle06 "[Tschumperle06]", 
/// and with background material in \ref anchor_Tschumperle02 "[Tschumperle02]".


//==============REFERENCES TO OTHER PAPERS AND RESOURCES=====================
///
/// \page sample_programs_page List of selected sample programs
///
/// \li <b>visualize the \ref anchor_StructureTensor "structure tensor" of an image</b>
///    <I>(compile and run \ref structureTensorComputation.C "app/structureTensorComputation.C") </I>\n
///    This small program will take an image as input, compute its structure tensor, and display the
///    result on screen.
///
/// \li <b>visualize application of the various filters used to discern image textures in segmentation</b>
///    <I>(compile and run \ref filter2DTests.C "app/filter2DTests.C") </I> \n
///    This small program will take an image as input, convert it to greyscale, and compute and visualize some filtered
///    versions of it that are useful for texture discrimination when segmenting images.  These are: 
///    \ref anchor_StructureTensor "structure tensor", anisotropically smoothed image (total variation variant)
///    and \ref anchor_ScaleMeasure "scale measure".
///
/// \li <b>generate a stack of 2D images that the library can read as a 3D volume</b> 
///    <I>(compile and run \ref generate_stack.C "app/generate_stack.C") </I>\n
///    This small utility program takes a list of image filenames, read the corresponding image files and generates
///    a "stack" of images that it writes to file in its own format.  A stack can be \ref lsseg::Image::read() "read" into
///    an \ref lsseg::Image "Image" as a 3D image, where each of the images in the stack constitutes one "z-slice" of 
///    the 3D volume.
///
/// \li <b>see the development of an image as TV-flow (total-variation-diminishing flow) is applied to it </b>
///     <I>(compile and run \ref tvfilter.C "app/tvfilter.C")</I> \n
///     This program will demonstrate the effect of total-variation-diminishing flow on an image.  Note the 
///     segmentation-like effect obtained.
///
/// \li <b>see the effect of \ref anchor_LIC "line integral convolution" on a 2D image</b> \n
///     <I>(compile and run \ref LIC_test.C "app/LIC_test.C")</I>\n
///     This small program demonstrates the effect of applying line integral convolution of a vector field and a 
///     Gaussian kernel function on an image consisting of random noise.  Note how this leads to an intuitive 
///     visualization of the vector field.
///
/// \li <b> see the effect of \ref anchor_LIC "line integral convolution" on a 3D image</b> \n
///     <I>(compile and run \ref field_smooth.C "app/field_smooth.C")</I>\n
///     This small program demonstrates the effect of applying line integral convolution of a Gaussian kernel function
///     and a 3D vector field generated by a dipole on a 3D volume consisting of random noise.  Note how this leads to an
///     intuitive visualization of the vector field.
///
/// \li <b>random partitioning of a rectangular domain into a given number of fragmented regions</b>
///     <I>(compile and run \ref randomVoronoiTest.C "app/randomVoronoiTest.C") </I> \n
///     Small program demonstrating how the function \ref lsseg::random_scattered_voronoi "random_scattered_voronoi()" can
///     be used to generate a random partitioning of a region into a given number of fragmented subregions.  This is
///     potentially useful for initializing an automatic segmentation process with multiple regions.
///
/// \li <b>2D segmentation with 2 or more regions regions</b> 
///     <I> (compile and run \ref segmentation2D.C "app/segmentation2D.C") </I>
///     This program carries out segmentation on a user-specified image according to various options that the user
///     can define on the command line.  Run the program without arguments to get a list of arguments to use.
///     The result from a segmentation can be saved afterwards.  Masks can also be used, defining which part(s) of the
///     image that should participate in the segmentation process.
///
/// \li <b>anisotropic, curvature-preserving smoothing of images</b> 
///    <I> (compile and run \ref smoothcurvpres.C "app/smoothcurvpres.C")</I>\n
///    This is a 2D/3D implementation of the \ref section_GreycStoration "greycstoration"
///    algorithm described by David Tschumperle in his 2006 paper
///    \ref anchor_Tschumperle06 "''Fast Anisotropic Smoothing of Multi-Valued Images using Curvature-Preserving PDE's''".
///    In essence, it is a curvature-preserving anisotropic smoothing algorithm using 
///    \ref anchor_LIC "line integral convolutions" to smooth an image along natural flow lines (edges, etc.)    
///    The results can be quite remarkable.  Try to play around with the parameters and see what it does.  An explanation
///    of the parameters is found in the paper linked to above.  The algorithm takes as input a picture in any common
///    image format, or a \em stack of images, which can be generated by the utility program
///    \ref generate_stack.C "generate_stack".
///
///
/// \page reference_page List of mentioned references
///
/// \anchor anchor_Brox05
/// <b>[Brox05]</b>\n
/// Thomas Brox: <I>From pixels to regions: partial differential equations in image analysis</I> \n
/// Faculty of Mathematics and Computer Science, Saarland University, Germany, april 2005 \n
/// Web link: 
/// <a href="http://www-cvpr.iai.uni-bonn.de/pub/pub/brox_PhDThesis.pdf">
/// http://www-cvpr.iai.uni-bonn.de/pub/pub/brox_PhDThesis.pdf</a> \n
///
/// \anchor anchor_Chan01
/// <b>[Chan01]</b>\n
/// T. F. Chan and L. A. Vese. <I>A level set algorithm for minimizing the Mumford-Shah
/// functional in image processing</I>  \n
/// In <I> Proc. First IEEE Workshop on Variational and Level Set Methods in Computer 
/// Vision </I>, pages 161-168, Vancouver, Canada, July 2001.  IEEE Computer Society
/// Press.\n
/// Web link: 
/// <a href="ftp://ftp.math.ucla.edu/pub/camreport/cam00-13.ps.gz">
/// ftp://ftp.math.ucla.edu/pub/camreport/cam00-13.ps.gz</a> \n
/// 
/// \anchor anchor_Chan99
/// <b>[Chan99]</b>\n
/// T. F. Chan and L. A. Vese: <I>An active contour model without edges</I> \n
/// In M. Nielsen, P. Johansen, O. F. Olsen and J. Weickert, editors, 
/// <I>Scale-Space Theories in Computer Vision</I>, volume 1682 of <I>Lecture Notes
/// in Computer Science</I>, pages 141-151. Springer, 1999. \n
/// Web link:
/// <a href="http://www.math.ucla.edu/~lvese/PAPERS/SS99.pdf">
/// http://www.math.ucla.edu/~lvese/PAPERS/SS99.pdf</a>\n
///
/// \anchor anchor_Cabral93
/// <b>[Cabral93]</b>\n
/// B. Cabral and L.C. Leedom.  <I>Imaging vector fields using line integral convolution</I>\n
/// SIGGRAPH'93, in Computer Graphics Vol.27, No.3, pp.263-272,1993
/// Web link:
/// <a href="http://portal.acm.org/citation.cfm?id=166151">http://portal.acm.org/citation.cfm?id=166151</a>
///  (subscription needed)\n
/// or <a href="http://www.cs.umu.se/kurser/TDBD13/VT00/extra/p263-cabral.pdf">
/// http://www.cs.umu.se/kurser/TDBD13/VT00/extra/p263-cabral.pdf</a>
/// 
/// \anchor anchor_Osher03
/// <b>[Osher03]</b>\n
/// S. Osher and R. Fedkiw: <I>Level set methods and dynamic implicit surfaces</I>\n
/// Volume 153 of Applied Mathematical Sciences.  Springer-Verlag, New York, 2003.\n
///
/// \anchor anchor_Tschumperle06
/// <b>[Tschumperle06] </b> \n
/// David Tschumperl&eacute;:
/// <I> Fast Anisotropic Smoothing of Multi-Valued Images using Curvature-Preserving PDE's </I> \n
/// International Journal of Computer Vision, IJCV(68), No 1, June 2006, pp. 65-82 \n
/// Web link:
/// <a href="http://www.greyc.ensicaen.fr/~dtschump/data/ijcv2006.pdf">
/// http://www.greyc.ensicaen.fr/~dtschump/data/ijcv2006.pdf</a>\n
///
/// \anchor anchor_Tschumperle02
/// <b>[Tschumperle02]</b> \n
/// David Tschumperl&eacute;:
/// <I>PDE-Based Regularization of Multivalued Images and Applications </I> \n
/// University of Nice-Sophia, December 2002 \n
/// FTP link:
/// <a href="ftp://ftp-sop.inria.fr/odyssee/Publications/PhDs/tschumperle:02.pdf">
/// ftp://ftp-sop.inria.fr/odyssee/Publications/PhDs/tschumperle:02.pdf</a>\n


//===================GLOSSARY OF SELECTED TERMS==============================
/// \page limited_glossary A glossary of a few selected terms as used in this documentation
///
/// \anchor anchor_ASmooth
///
/// \li <b> Anisotropic smoothing</b> \n
/// In short, anisotropic smoothing of an image, as opposed to \em isotropic smoothing,
/// is some kind of regularization process that does not smooth the image equally in all
/// directions.  One of the earliest schemes for anisotropic image smoothing was 
/// Perona and Malik's regularization, which you can read a short presentation about
/// in section 2.1.2 of \ref anchor_Tschumperle02 "D. Tschumperle's PhD thesis".
/// 
/// \anchor anchor_Channel
/// \li By an image <b>channel</b> \n
/// we mean one of the image components.  Each pixel in the image contains one
/// separate value for each image channel.  For instance, a color image usually
/// has at least three separate channels, one for its red component, another for its
/// green and a third for its blue.  Channels do not conceptually need to represent
/// colors though.  the \ref lsseg::LevelSetFunction "LevelSetFunction" class, which
/// inherits from the \ref lsseg::Image "Image" class, always contains exactly \em one
/// channel, which represents the value of the \ref anchor_LevelSetFunction "level-set function"
///  for each pixel position.
/// A \ref lsseg::Mask "Mask" is another kind of single-channeled image, where the
/// channel value per-pixel is boolean and represents whether or not that pixel is \em active
/// or <em> masked out </em>.
///
/// \anchor anchor_DiffusionTensor
/// \li A <b> diffusion tensor </b> \n
/// in 2D is a tensor described by a symmetric, semi-definite, positive 2x2 matrix (3x3 in 3D).
/// It is used in smoothing partial differential equations, and describes the amount of diffusion
/// along privileged spatial directions.  The matrix has two real eigenvalue/eigenvector pairs,
/// its eigenvectors are perpendicular.  The eigenvectors describe two spatial directions, and their
/// respective eigenvalues represent the amount of diffusion along these directions.
/// If we denote its eigenvalues \f$\lambda_1, \lambda_2\f$ and the corresponding eigenvectors
/// \f$u_1, u_2\f$, the diffusion tensor \f$T\f$ can be decomposed into:
/// \f[T = \sum_{i=1}^2 \lambda_i u_i u_i^T\f]
///
/// \anchor anchor_EulerLagrange
///
/// \li The <b> Euler-Lagrange equation</b> \n
/// The Euler-Lagrange equation is the major formula of the calculus of variations.  It can
/// be considered the gradient of a functional, and gives a necessary condition that must be 
/// verified for a function \f$f\f$ that minimizes the functional.
/// For functions over an \f$n\f$-dimensional domain \f$\Omega\f$, and for a functional 
/// \f[J(f)=\int_\Omega F(f, x_1, ..., x_n, f_{x_1}, ..., f_{x_n}) d\Omega\f]
/// the corresponding Euler-Lagrange equation is given by:
/// \f[\frac{\partial F}{\partial f} - \sum_{i=1}^n \frac{\partial}{\partial x_i}\frac{\partial F}{\partial f_{x_i}} = 0\f]
/// As the right-hand side of the equation above can be considered the gradient of \f$J\f$, 
/// we can use it in a way analog to \em function minimization by gradient descent.  In this
/// way, we look for a minimum of \f$J\f$ by iteratively taking steps going "downhill" until the
/// gradient vanishes.  We introduce an "artificial" time parameter, then start from an initial function \f$f_0\f$ and
/// following the opposite direction of the gradient until we arrive at a local minimum of \f$J\f$.  
/// This leads to the following partial differential equation evolution process:
/// \f{eqnarray*}
        f_{t=0} &=& f_0 \\ 
	\frac{\partial f}{\partial t} &=& - \left(\frac{\partial F}{\partial f} - \sum_{i=1}^n\frac{d}{d x_i}\frac{\partial F}{\partial f_{x_i}}\right) \\
   \f}
///
/// \anchor anchor_Functional
///
/// \li A <b>functional</b> \n
/// is, for our purposes, a function that takes functions (from some specified class) as its argument, and return a
/// scalar value.  We are often interested in the functions for which a  given functional reaches a minimum.  
/// It is typically defined by an integral.  For functions over a domain \f$\Omega\f$, this can be written
/// \f[J(f) = \int_\Omega F(x, f(x), \nabla f(x), ...)d\Omega\f] 
/// where \f$F(...)\f$ is a function of the position in \f$\Omega\f$, the value of \f$f\f$ and its derivatives
/// (possibly also higher-order) at this position.
/// In the process of searching for a minimum of a given functional, we often make use of the
/// \ref anchor_EulerLagrange "Euler-Lagrange equation".
///
/// \anchor anchor_HeatEquation
///
/// \li The <b> heat-equation </b> \n
/// is a partial differential equation that is the simplest form of the diffusion equation.  
/// It is much used in physics to describe diffusion of a quantity through an isotropic medium, like \em heat
/// through a \em solid.  If \f$I\f$ is the function describing the quantity, the equation can be written:
/// \f{eqnarray*}  I_{(t=0)} & = & I_0 \\ \frac{\partial I}{\partial t} & = & \Delta I \f}
///
/// \anchor anchor_Hessian
/// 
/// \li the <b> Hessian matrix \f$H\f$</b> \n
/// of a function \f$f : R^2 \rightarrow R\f$ is defined by the second derivatives and 
/// cross derivatives of \f$f\f$:
/// \f[H(x,y) := \left[ \begin{array}{cc} \frac{\partial^2 f}{\partial x^2} & \frac{\partial^2 f}{\partial x \partial y} \\ \frac{\partial^2 f}{\partial x \partial y} & \frac{\partial^2 f}{\partial y^2} \end{array} \right]  \f]
/// \anchor anchor_LevelSetFunction
/// 
/// \li A <b>level-set function</b> \n
/// is, in our terminology, a scalar function \f$\phi\f$ defined on some domain \f$\Omega\f$ and which is used to define
/// a partitioning of the domain into two regions separated by a boundary defined by its \ref anchor_ZeroSet "zero set".
/// Hence, \f$\phi\f$ takes on positive value for one region, and negative for the other.  The boundary curve is
/// \em implicitly \em defined by the level-set curve.  On the other hand, for any given closed curve in \f$\Omega\f$, we
/// can define a level-set function that represents it.  A common way of doing that is by letting the level-set function
/// be the \ref anchor_SignedDistanceFunction "signed distance function" of the curve. Level-set functions is defined
/// in LSSEG by the class \ref lsseg::LevelSetFunction "LevelSetFunction", which is a very central data structure.
///
/// \anchor anchor_LIC
///
/// \li <b> Line-Integral-Convolution (LIC) </b> \n
/// Given an image, a vector field defined on the same domain \f$\Omega\f$ as the image, and 
/// a kernel function \f$k : \mathbf{R} \rightarrow \mathbf{R}\f$, we define the line-integral-convolution
/// as the one-dimensional convolution of the image with the kernel function along all streamlines
/// defined by the vector field.  This results in an image that is smoothed along those streamlines.
/// This technique is useful for visualizing a vector field (by convoluting an image consisting only of noise
/// with a gaussian kernel along the streamlines of the field), but also for curve-preserving smoothing
/// schemes of images (in which case the vector field has been determined based on structures present
/// in the original image).  Refer to \ref anchor_Cabral93 "[Cabral93]" for details on the process, and to
/// \ref anchor_Tschumperle06 "[Tschumperle06]" for this technique used for curve-preserving smoothing.
///
/// \anchor anchor_MeanCurvatureMotion
///
/// \li <b> Mean curvature motion </b> \n
/// is a motion of an interface (boundary curve, surface) along its normal direction, with a velocity proportional to its
/// curvature.  It can generally be described by the equation:
/// \f[\overrightarrow{V} = -b\kappa \overrightarrow{N}\f]
/// where \f$\overrightarrow{V}\f$ is the velocity vector, \f$b > 0\f$, \f$\overrightarrow{N}\f$ is the normal vector
/// and \f$\kappa\f$ is the curvature.
/// In a level-set formulation, where the interface is described by a \ref anchor_LevelSetFunction "level-set function",
/// this equation can be written:
/// \f[\frac{\partial\phi}{\partial t} = b\kappa |\nabla\phi|\f]
/// where \f$\phi\f$ is the level-set function and \f$\kappa = \nabla . \left(\frac{\nabla\phi}{|\nabla\phi|}\right)\f$ is the 
/// curvature of the isocurves of \f$\phi\f$.  The full formulation is therefore:
/// \f[\frac{\partial\phi}{\partial t} = b \nabla . \left(\frac{\nabla\phi}{|\nabla\phi|}\right) |\nabla\phi|\f]
/// Numerically, this can be discretized using <em> central differencing </em> (as decribed by \ref anchor_Osher03 "[Osher03]".
/// In order for the integration to be stable, the timestep \f$\Delta t\f$ must respect the following CFL-condition (on a 2D domain):
/// \f[\Delta t\left(\frac{2b}{(\Delta x)^2} + \frac{2b}{(\Delta y)^2}\right) < 1\f]
/// Motion by mean curvature is described in chapter 4 of \ref anchor_Osher03 "[Osher03]".
///
/// \anchor anchor_MumfordShahFunctional
///
/// \li The <b> Mumford-Shah functional </b> \n
/// is an energy \ref anchor_Functional "functional" used in the context of image 
/// segmentation.  For a given partitioning of an image into region, and
/// a given approximation of the original image in those regions, it 
/// gives a numerical value that measures the "goodness" of that particular
/// partitioning.  This value takes into account the difference between
/// the approximated image and the original one, as well as how regular
/// the proposed approximation is and how long the combined length of boundaries.  
/// Lower values of this functional means a better segmentation, and for that reason, 
/// image segmentation often becomes the matter of searching for minima of this functional.  
/// In its full form, it can be written:
/// \f[E(u, \Gamma) = \int_\Omega (u-I)^2 dx + \lambda \int_{\Omega - \Gamma} |\nabla u|^2 dx + \nu\int_\Gamma ds\f]
/// Here, \f$\Gamma\f$ represent the partitioning of the image domain \f$\Omega\f$ into distinct regions,
/// whereas \f$u\f$ is a smooth approximation of the image within each region (it is allowed to be discontinuous 
/// across region boundaries, though).  The first term on the right side measures the distance between the approximated
/// image \f$u\f$ and the origional image \f$I\f$.  The second term on the right side measures the regularity of \f$u\f$ 
/// within each region.  The last term on the right side measures the total length of the boundaries separating the 
/// regions. \f$\lambda\f$ and \f$\nu\f$ are tuneable weighing terms describing how much importance to give to
/// the second and third right-hand-side term compared to the first.
///
/// \anchor anchor_NormalDirectionFlow 
///
/// \li The <b> normal direction flow </b> \n
/// is referred to here as a motion of an interface (boundary curve, surface) along its normal direction, vith a velocity 
/// imposed from outside.  We suppose that this motion does not depend on the interface, but it \em can depend on 
/// the spatial coordinate.  
/// In the level-set formulation, where the interface is described by the zero-set of a
/// \ref anchor_LevelSetFunction "level-set function" \f$\phi\f$, the formula is:
/// \f[\frac{\partial\phi}{\partial t} + a|\nabla\phi| = 0\f]
/// (Here, \f$a\f$ might vary over the domain, although it is not dependent on \f$\phi\f$).
/// This is an example of a <em> Hamilton-Jacobi </em> equation, and must be discretized accordingly.
/// In order to assure stability, the timestep used must then fulfill the following CFL condition (on a 2D domain):
/// \f[\Delta t \max\left\{\frac{|H_1|}{\Delta x} + \frac{|H_2|}{\Delta y}\right\} < 1\f]
/// where \f$H_1\f$ and \f$H_2\f$ are the partial derivatives of the systems Hamiltonian \f$H\f$ with respect to 
/// \f$\phi_x\f$ and \f$\phi_y\f$.  (The Hamiltonian of this equation is \f$ H(\nabla\phi) = a|\nabla\phi|\f$).
/// Motion in the normal direction is described in chapter 6 of \ref anchor_Osher03 "[Osher03]", although there, \f$a\f$
/// does not vary over the domain.  The theory of Hamilton-Jacobi equations and their numerical discretization
/// is introduced in chapter 5 of the same book.
///
///
///  \anchor anchor_NormalForceField
/// \li
/// A <b> normal force field </b> is a field of scalars defined over the domain of a
/// level-set function \f$\phi\f$, representing a force acting normally on the level-set
/// curve (surface in 3D) of \f$\phi\f$ at that particular point.  Normal force fields are 
/// used to evolve level-set curves (surfaces) over time, and are determined using
///  \ref lsseg::ForceGenerator "ForceGenerators". 
/// \note Since a normal force field can be seen as a single-valued image over the domain, 
/// in the LSSEG library we represent it using another \ref lsseg::LevelSetFunction 
/// "LevelSetFunction" (\em not to be confused  with the \ref lsseg::LevelSetFunction
/// "LevelSetFunction" used to represent the actual segmentation).
///
/// \anchor anchor_Reinitialization
/// 
/// \li <b>Reinitialization</b> of a level-set function \n
/// During the segmentation process, the level-set function undergoes an evolution based 
/// on a partial-differential equation.  What interests us is where the level-set function
/// is negative, where it is positive and where it is zero (i.e. its
/// \ref anchor_ZeroSet "zero set").  Beyond that, its exact numerical values are generally
/// not of interest for us.  However, during the evolution, the level-set function can 
/// develop regions that are very flat or very steep, which might stall the process due to
/// increasingly shorter stable timesteps, or other numerical problems.  For that reason
/// it is desireable to reset the level-set function to a
/// \ref anchor_SignedDistanceFunction "signed distance function" once in a while.  The
/// signed distance function in question should have same zero-set and sign as the original
/// level-set function, but having a gradient of 1 almost everywhere, it leads to a much more 
/// stable evolution process afterwards.
///
/// \anchor anchor_ScaleMeasure
///
/// \li The <b> scale measure </b> \n
/// as defined by Thomas Brox in section 4.1.2 of \ref anchor_Brox05 "his thesis", is used
/// in the context of texture discrimination for image segmenting purposes.  It is a 
/// single-channeled image derived from the image to be segmented, and whose pixel values give
/// an estimate of the \em scale of the structure that the corresponding pixels in the original
/// image are part of.  The computation of this measure is based on a 
/// total-variation-diminishing flow, and is somewhat technical.  The reader is referred to 
/// the relevant section of \ref anchor_Brox05 "Brox' thesis" for all the details.
///
/// \anchor anchor_SignedDistanceFunction
/// \li
/// The <b> signed distance function </b> of a closed contour in a 2D domain
/// or a closed surface in a 3D domain, is the function that assigns to each 
/// point in the domain <em>its shortest distance to the contour/surface</em>, with 
/// a sign that indicates whether the point is situated \em inside or \em outside
/// the enclosed domain.  In our implementation, we have chosen the convention that
/// negative values signifies the \em inside.  The \ref anchor_ZeroSet "zero-set"
/// of a signed distance function is exactly the closed contour/surface we started
/// out with.
///
/// \anchor anchor_StructureTensor
/// \li The <b> structure tensor</b> \n
/// of an image \f$I\f$ is a tensor field over the image domain \f$\Omega\f$.
/// In each point of \f$\Omega\f$, the tensor describes the orientation of local
/// structures in the image around this point.  It is represented by a symmetric,
/// positive, semi-definite 2x2 matrix.  For \em greyscale images, this matrix, \f$G_0\f$, is
/// defined by (subscripts of \f$I\f$ denote partial derivatives):
/// \f[G_0(x,y) := \nabla I\nabla I^T = \left[ \begin{array}{cc} I_x^2 & I_xI_y \\ I_xI_y & I_y^2 \end{array} \right]  \f]
/// For multi-channel images (e.g., color images) with N channels, the definition is:
/// \f[G_0(x, y) := \sum_{i=1}^N \nabla I_i \nabla I_i^T \f]
/// For various reasons, we often prefer to work with the \em smoothed structure tensor \f$G_\rho\f$ instead.
/// This is obtained through convolution with a Gaussian kernel:
/// \f[G_\rho : = G_0 * K_\rho\f]
/// Much useful information about the structure tensor can be found in section 2.2 of
/// \ref anchor_Brox05 "Brox's PhD thesis".
/// 
///
/// \anchor anchor_ZeroSet
/// \li
/// The <b> zero-set </b> of a level-set function is the ensemble of points in 
/// its domain for which the function evaluates to zero.


//=================== PAGE ON HOW LSSEG SEGMENTATION WORKS ==============================
///
/// \page page_HowSegmentationWorks How Segmentation in LSSEG Works
///
/// The theory behind the level-set based segmentation approach used by LSSEG is thoroughly
/// explained in chapter 4 and 5 of \ref anchor_Brox05 "Thomas Brox' PhD thesis".  
/// The reader is strongly recommended to read at least these two chapters, as the passages
/// below only are intended to give some general idea about what is going on.
///
/// \section section_OverallSegmentationExplanation Overall explanation of problem
///
/// When we talk about \em segmentation of an image, we mean a way to partition the image
/// into different regions in some meaningful way.  Generally, we want to separate parts
/// of the image that belong to different objects.   When attempting to do this in an
/// automatical fashion, we run into the problem that a computer does not have any direct
/// way of recognizing physical objects.  For that reason, we define segmentation algorithms
/// that divide the image into regions which are similar in some specifiable way, and try
/// to make this correspond as well as possible with the actual image semantics.
///
/// The way we work, we try to define regions whose boundaries are somewhat regular, and 
/// in a way such that the pixel distribution of each region has distinct statistical 
/// properties.
///
/// \section section_LevelSetRepresentation How regions are represented - level set functions
/// 
/// If we want to separate an image into regions, we need some way to represent a region, ie.
/// a way to specify which of the image pixels that are contained within the region.  
/// A region does not necessarily need to be connex (it can consist of multiple disjoint parts)
/// and its topology can include handles and holes.
/// 
/// In this library, we use \ref anchor_LevelSetFunction "level-set functions" to represent
/// regions.  They are very flexible in that they allow for arbitrary topology of regions,
/// are easy to modify and as (differentiable) functions they fit well into a PDE framework.
///
/// A level-set function is a scalar function defined over the whole image domain, and defines a 
/// separation of the image into two parts: those parts of the image where the level-set function
/// takes on a negative value is defined as being a separate region from those parts where the
/// level-set function takes on a positive value.  In LSSEG, we choose the convention that all
/// parts of the image domain where the level-set function is \em negative is called 
/// the \em interior of the region represented by the level-set function (all the rest thus
/// constituting the \em exterior).
/// Also, those parts of the image where the level-set function is \em zero is called
/// the \em boundary of the domain.
///
/// in LSSEG, the level-set function is represented by the class
///  \ref lsseg::LevelSetFunction "LevelSetFunction".  This class inherits from the
///  \ref lsseg::Image "Image" class and can thus be regarded as a (greyscale) image itself;
/// each of its pixels representing the value taken on by the level-set function at the position
/// of that pixel.  (NB: the \ref lsseg::LevelSetFunction "LevelSetFunction" has a restriction
/// not present in its base class \ref lsseg::Image "Image": it always has exactly \em one
/// \ref anchor_Channel "channel").
///
/// For a very good book on level-set functions and their use, refer to
/// \ref anchor_Osher03 "[Osher03]".
/// \subsection subsection_ChangingRegionShape changing the shape of a region described by a level-set function
///
/// In the segmentation process described here, we obtain our final partitioning of the image by starting
/// from some initial partitioning, which we change in a continuous way in order to arrive at the region(s)
/// (hopefully) sought.  This is done by developing the level-set function describing a region, 
/// according to a time-dependent partial differential equation that incorporates information taken from 
/// the image.
/// So how can we change a region using a time-dependent partial differential equation?  We will here describe
/// two basic equations we use, and the effect that they have on a region described by a level-set function \f$\Phi\f$.
/// Remember our convention that the interior of the region described by our level-set function is the domain
/// on which the level-set function is \em negative, and that the \em boundary of this region is described by
/// the domain on which the level-set function is exactly zero (a curve, except in degenerate cases).
/// Again, refer to the \ref anchor_Osher03 "book by Osher" for formal explanations and additional detail.
///
/// \li <b>mean curvature motion</b>
/// The main effect of mean curvature motion is to shrink the interior region and smoothen its boundary curve.
/// If mean curvature alone is allowed to act on a (closed) region described by a level-set function, this region
/// will eventually shrink until it disappears.  The partial-differential equation is: 
/// \f$\partial_t \Phi = |\nabla\Phi|div(\frac{\nabla\Phi}{|\nabla\Phi|})\f$.  Discretized, this describes the
/// change in \f$\Phi\f$ from one timestep to the next.  In the LSSEG library, the quantity
/// \f$ |\nabla\Phi|div(\frac{\nabla\Phi}{|\nabla\Phi|}) \f$ can be directly obtained from a
/// \ref lsseg::LevelSetFunction "LevelSetFunction" through use of the member function 
/// \ref lsseg::LevelSetFunction::curvatureTimesGrad2D() "LevelSetFunction::curvatureTimesGradXD()".
/// Note that movement by mean curvature motion is completely determined by the shape of the level-set function
/// itself, and is not influenced by any external factors.
///
/// \li <b>\ref anchor_NormalDirectionFlow "normal direction flow" </b>
/// The effect of normal direction flow is to push the boundary curve outwards / pull the boundary curve
/// inwards along its normal in each point.  The magnitude of this "push" can vary over the domain, and 
/// in our segmentation process, this magnitude is governed by information taken from the picture.
/// The partial differential equation is: \f$\partial_t\Phi = \beta \dot |\nabla\Phi|\f$
/// Here, \f$\beta\f$ may vary over the domain, hence the effect on \f$\Phi\f$ is not the same everywhere,
/// which allows for controlling its development according to external factors.
/// The incremental change in the level-set function due to normal direction flow can be computed in 
/// the LSSEG library by using the functions \ref lsseg::normal_direction_flow_2D() "normal_direction_flow_XD()".
/// In order to define a \f$\beta(x,y)\f$ function to use in segmentation, we use \ref lsseg::ForceGenerator "ForceGenerators".
///
/// <span style="color: red">In the segmentation process of this library, these are the two only kinds of
/// movement tools necessary. </span>
///
/// \note \li It might be not be intuitive to associate the movement of a boundary curve with the development
/// of a level-set function.  In order to see this clearer, the reader should consult the 
/// \ref anchor_Osher03 "book by Osher", or read section 5.1 of \ref anchor_Brox05 "Brox' thesis".
///
/// \note \li Applying any of the above differential equations on \f$\Phi\f$ is quite computationally expensive.
/// However, we are only interested in the development of \f$\Phi\f$ around its current zero-set, where the 
/// current region boundary is located.  Therefore, we only apply the above equations to a small band along
/// this boundary.  Since this boundary will move over time, the band will move along with it.  All parts
/// of the level-set function away from this band will be left as it is (except that it may be occationally
/// \ref anchor_Reinitialization "reinitialized" to the \ref anchor_SignedDistanceFunction "signed distance function").
/// 
/// 
///
/// \section section_MumfordShah Specifying a segmentation criterion - the Mumford-Shah functional
///
/// In order to partition our image into regions, we need to specify what we are trying to obtain.
/// As we do not assume any pre-knowledge of the various objects that a picture can contain, we must
/// establish some more general criteria.  For image segmentation, is usual to depart from the
/// purely mathematical construct known as the <em> Mumford-Shah functional </em>.  It looks like this:
///
/// \f[ E(u, \Gamma) = \int_\Omega (u-I)^2 dx + \lambda \int_{\Omega - \Gamma}|\nabla u|^2 dx + \nu \int_\Gamma ds\f]
///
/// In the above formula, \f$\Omega\f$ represents the whole image domain, whereas \f$\Gamma\f$ represents
/// the boundaries that make up an image partition.  \f$u\f$ is a simpler image, which we can consider as 
/// "a way to approximate the contents of each image region".  
///
/// The Mumford-Shah is a \ref anchor_Functional "functional" that assigns one numerical value for each 
/// possible partition \f$\Gamma\f$ of the image domain (with a given way \f$u\f$ to approximate the contents 
/// of each region), and it is chosen in such a way that "better" partitions are assigned lower values.
/// ("Better" in the sense that the bondaries are kept short and that each domain has distinct properties that
/// can be modeled well by out simpler approximation \f$u\f$).  Thus, in order to search for a good
/// segmentation, we search for a partitioning of the image that gives a value that is as low as possible
/// when applying the Mumford-Shah functional.  You can 
/// \ref anchor_MumfordShahFunctional "read more about the Mumford-Shah functional in the glossary".  
/// In order to get a more formal and in-depth explanation about how it is used for our segmentation purposes, 
/// you should consult chapter 5 of \ref anchor_Brox05 "[Brox05]".
///
/// The Mumford-Shah functional provides a very general setting for describing image segmentation.  Many 
/// segmentation methods can be seen as different ways of searching for minima of this functional, with different
/// ways of representing \f$\Gamma\f$ and different approximations \f$u\f$.  In our PDE-based approach, we are 
/// interested in the functional's derivatives, which we suppose exist, and to use them in a minimization process 
/// in order to converge to a minimum, which we define as \em the segmentation we are looking for.
///
/// \section section_ObtainingPDE Obtaining our partial differential equation
/// 
/// The gradient of the above functional is determined by the \ref anchor_EulerLagrange "Euler-Lagrange equations".
/// The exact form of the obtained equations depend on how we represent \f$\Gamma\f$ (region boundaries) and \f$u\f$ 
/// (approximating model of each region).  In our setting, for now, we aim to partition the image into <em>exactly two</em>,
/// <em>not necessarily connex</em>, regions, and represent this partitioning by a level-set function \f$\Phi\f$.  
/// \f$\Gamma\f$ is then the zero-set of \f$\Phi\f$.
///
/// \li <b>(Chan-Vese model)</b>: We can define \f$u\f$ to be piecewise constant on each of the two regions defined by
///  \f$\Phi\f$.  This means that
/// the second term of the Mumford-Shah functional above vanish (\f$\nabla u\f$ zero everywhere), and that the constant
/// value of \f$u\f$ for a given region must be equal to the average pixel value for that region, in order to minimize the
/// functional's first term.  Through the Euler-Lagrange equations we obtain an expression for an evolving \f$\Phi\f$ that, 
/// in a slightly modified version, can be written:
/// \f[\partial_t\Phi = |\nabla\Phi|\left[\nu\nabla(\frac{\nabla\Phi}{|\nabla\Phi|})-(I - u_{outside})^2 + (I - u_{inside})^2\right]\f]
/// \f$I\f$ here represents the image, so it varies over the domain.  \f$u_{outside}\f$ and \f$u_{inside}\f$ are the constant
/// values for \f$u\f$ inside and outside the specified region.  These values are obtained by averaging pixel values.
/// This is the segmentation model proposed by \ref anchor_Chan99 "Chan and Vese in 1999".  For details on how we 
/// arrive at this equation, you can also consult \ref anchor_Osher03 "chapter 12 of Osher's book", or subsection 5.1.2 of
/// \ref anchor_Brox05 "Brox' PhD thesis".  
/// We note that, by posing \f$\beta = (I - u_{outside})^2 + (I - u_{inside})^2\f$, the above expression can be rewritten:
/// \f[\partial_t\Phi = \nu |\nabla\Phi| \left[\nabla(\frac{\nabla\Phi}{|\nabla\Phi|})\right] + \beta |\nabla\Phi|\f]
/// As you can see the right-hand side now consists of two terms which are the expressions for
/// \ref anchor_MeanCurvatureMotion "mean curvature motion" and \ref anchor_NormalDirectionFlow "normal direction flow"
/// respectively, with \f$\nu\f$ playing the role of deciding the weight of one term with respect to the other.
/// As stated earlier, in LSSEG we use \ref lsseg::ForceGenerator "ForceGenerators" in order to compute \f$\beta(x,y)\f$
/// for various segmentation schemes.  The ForceGenerator in LSSEG corresponding to the Chan-Vese model is called 
/// \ref lsseg::ChanVeseForce "ChanVeseForce".
///
/// \anchor anchor_RegionalStatistics
/// \li <b>(Regional statistics)</b>: More elaborate segmentation models can be obtained if we allow more complicatepd
/// expressions for the approximation \f$u\f$ in the \ref anchor_MumfordShahFunctional "Mumford-Shah functional".  n
/// If we let \f$u\f$ represent some probability distribution of pixels, we obtain the following equation of evolution
/// for \f$\Phi\f$:
/// \f[\partial_t\Phi = |\nabla\Phi|\left[\log p_{outside} - \log p_{inside}\right] + \nu|\nabla\Phi|\left[\nabla(\frac{\nabla\Phi}{|\nabla\Phi|})\right] \f]
/// Here, \f$p_{outside}\f$ and \f$p_{inside}\f$ represent, respectively, the probability of a particular pixel value inside
/// and outside of the region defined by our level-set function \f$\Phi\f$.  We see that, again, the equation of evolution
/// is a sum of a term describing \ref anchor_MeanCurvatureMotion "mean curvature motion" and a term describing
/// \ref anchor_NormalDirectionFlow "normal direction flow", with \f$\nu\f$ being a parameter that determines the relative
/// weight of the normal direction flow term with respect to the other term.  Read subsection 5.1.2 of
/// \ref anchor_Brox05 "Brox' PhD thesis" in order to see how we arrive at this expression. 
/// \n
/// If we use a normal probability distribution to model pixel distributions, the above 
/// expression becomes:
/// \f[\partial_t\Phi = |\nabla\Phi|\left[\frac{(I - \mu_{inside})^2}{2\sigma_{inside}^2} - \frac{(I - \mu_{outside})^2}{2\sigma_{outside}^2}\right] + \nu|\nabla\Phi|\left[\nabla(\frac{\nabla\Phi}{|\nabla\Phi|})\right] \f]
/// where \f$\mu_{inside}\f$ and \f$\mu_{outside}\f$ are the estimated average pixel values of the image inside and outside of
/// the region, and \f$\sigma_{inside}^2\f$ and \f$\sigma_{outside}^2\f$ are the estimated variances.  Contrary to the
/// <b>Chan-Vese model</b>, this model is able to distinguish between two regions having the same average pixel values, as long
/// as there is a distinct difference in the variances of the two regions.  The LSSEG \ref lsseg::ForceGenerator "ForceGenerator"
/// corresponding to this model is \ref lsseg::NormalDistributionForce "NormalDistributionForce".
/// We can choose other probability distributions as well; if we construct probability distributions based on the current 
/// histograms of pixels in each region, the Parzen density estimate, we obtain an even richer model, but which is more 
/// dependent on a good initalization.  The LSSEG \ref lsseg::ForceGenerator "ForceGenerator" used in this model is called
/// \ref lsseg::ParzenDistributionForce "ParzenDistributionForce".  Again, for details, refer to 
/// \ref anchor_Brox05 "Brox' PhD thesis".
///
/// The algorithm in LSSEG that carries out the evolution of a level-set function \f$\Phi\f$ based on a selected 
/// \ref lsseg::ForceGenerator "ForceGenerator", \ref lsseg::Image "Image", initialization and smoothness value \f$\nu\f$ 
/// is called \ref develop_single_region_2D() "develop_single_region_XD()".
///
/// \section section_TextureAndMultichannels Recognizing texture and working with multiple channels
///
/// In the setting above, we supposed that the image \f$I\f$ was greyscale, i.e., contained only one channel.  It is however
/// very easy to extend this approach to a multi-channel setting.  In the evolution equation, the term expressing 
/// \ref anchor_MeanCurvatureMotion "mean curve motion" remains the same, whereas the term expressing 
/// \ref anchor_NormalDirectionFlow "normal direction flow" is modified to take into account the contribution from each
/// channel   For the model using regional statistics, we will then obtain (supposing we have M different, independent
/// channels):
/// \f[\partial_t\Phi = |\nabla\Phi|\left[\sum_{j=1}^M(\log p_{outside, j} - \log p_{inside,j})\right] + \nu|\nabla\Phi|\left[\nabla(\frac{\nabla\Phi}{|\nabla\Phi|})\right] \f]
/// Details are found in section 5.2 of \ref anchor_Brox05 "Brox' PhD thesis".  
/// The LSSEG \ref lsseg::ForceGenerator "ForceGenerators" are designed to work with multichanneled images (including, of course
/// the special case of single-channeled images).
///
/// No matter which of the above mentioned models we use, we still need to ensure that the regions into which we want to
/// separate the image have clearly discernable statistical properties.  For a white mouse on a black carpet, this is 
/// easy - the average pixel intensity is very different between the "mouse" region and the "carpet" region, and a 
/// greyscale, single-channeled image should be sufficient to obtain a good result.  For an image describing a red cup on a
/// blue table, it is likewise a fairly easy task to obtain the segmentation by using (at least) the "RED" and "BLUE" color
/// components of the color picture.
/// However, often we would like to segment regions that are not easily distinguishable by color or luminosity.  In this
/// case, we would have to base our segmentation on different texture properties of the regions we want to obtain.  For this
/// purpose, we must find ways to filter the image in ways that enhance the pixel distribution statistical properties between
/// regions of different textures.  This is not an easy thing to do, but two possibilities are:
///
/// \li the \ref anchor_StructureTensor "structure tensor": (see link for definition)\n
/// This actually creates three different channels, based on the image's x-derivative, y-derivative and cross-derivative.
/// In LSSEG, these channels can be computed using the function
/// \ref compute_structure_tensor_2D() "compute_structure_tensor_XD()".
/// Click <a href="../images/zebraTensor.png" target="_blank">here</a> for a higher-resolution version of the below image.\n\n
/// \image html zebraTensor_small.png "An image of a zebra, with its three structure tensor components."
/// You can read more about the structure tensor in section 2.2 of \ref anchor_Brox05 "Brox' PhD thesis".
///
/// \li the \ref anchor_ScaleMeasure "scale measure": (see link for definition)\n
/// The scale measure attempts to express the size of the structures each pixel belong to, in the sense that pixels that have
/// the same scale measure value are likely to belong to structures with the same local scale.  
/// In LSSEG, the scale measure can be computed using the function
/// \ref compute_scale_factor_2D() "compute_scale_factor_XD()".
/// Click <a href="../images/zebraScale.png" target="_blank">here</a> for a higher-resolution version of the below image. \n\n
/// \image html zebraScale_small.png "An image of a zebra, with its computed scale measure (computed using 20 timesteps of duration 1)"
/// You can read more about the scale measure in section 4.1.2 of \ref anchor_Brox05 "Brox' PhD thesis".
///
///
/// \section section_MultiRegion Segmentation with multiple regions
///
/// The LSSEG library supports segmentation of a picture into more than one region.  This is a complex problem, and several
/// approaches to this exist in the literature.  We have implemented an algorithm similar to what is described in section
/// 5.3 of \ref anchor_Brox05 "Brox' PhD thesis".  Some characteristics of this approach are:
/// \li we use one separate level-set function for each region.  This is somewhat different from the usual two-region 
/// approach, where \em one level-set function was sufficient to determine both regions ("inside" and "outside").
/// \li the number of regions must be specified in advance by the user.  
/// \li there should be no substantial overlap between the segmented regions.  For algorithmic reason, an overlap of a few 
/// pixels is tolerated.
/// \li whereas some multiregion segmentation approaches found in the literature use <em>lagrange multipliers</em> in order
///  to keep different regions separate, we use a <em>competing region</em> approach, where the regions compete for influence
/// along their mutual borders, and can "push" each other away.
///
/// Whereas the theoretical development can be found in section 5.3 of \ref anchor_Brox05 "Brox' PhD thesis", we here aim
/// to outline the general idea and what consequences this leads to in the LSSEG setting.
///
/// In the original evolution equation for segmentation into two regions using
/// \ref anchor_RegionalStatistics "regional statistics", we had:
/// \f[\partial_t\Phi = |\nabla\Phi|\left[\log p_{outside} - \log p_{inside}\right] + \nu|\nabla\Phi| \left[\nabla(\frac{\nabla\Phi}{|\nabla\Phi|})\right] \f]
/// In the multi-region setting, there is no \em outside anymore, just a number \f$N\f$ of different competing regions described
/// by \f$N\f$ level-set functions \f$\Phi_i, i=1..N\f$.
/// For a given level-set function \f$\Phi_i\f$, we modify its evolution equation in the following way:
/// \f[\partial_t\Phi_i = |\nabla\Phi_i|\left[\log p_i + \nu div(\frac{\nabla\Phi_i}{|\nabla|Phi_i|}) - \max_{j\neq i} (\log p_j + \nu div(\frac{\nabla\Phi_j}{|\nabla\Phi_j}))\right]\f]
///
/// We see that the term representing the outside (\f$ \log p_{outside}\f$) in the original equation has now been replaced by
/// a term \f$\log p_j + \nu div(\frac{\nabla\Phi_j}{|\nabla\Phi_j})\f$ representing the "strongest" competing candidate at
/// the point in question.  This leads to an evolution where stronger candidates "push away" weaker ones.  As long as we take
/// care to initialize the domains so that they are non-overlapping and together cover the whole image domain, this will lead
/// to an evolution where regions stay separate while still covering the whole image domain at all times.  
///
/// The multi-region setting has one important implication for the initialization of \ref lsseg::ForceGenerator "ForceGenerators"
/// when using the LSSEG library: by default the ForceGenerators \ref lsseg::NormalDistributionForce "NormalDistributionForce"
/// and \ref lsseg::ParzenDistributionForce "ParzenDistributionForce" compute the quantity corresponding to the
/// \ref anchor_NormalDirectionFlow "normal direction flow" for \em two-region segmentation, i.e., 
/// \f$(\log p_{outside} - \log p_{inside})\f$.  In the multi-region setting, the \em outside is not defined anymore, and the
/// corresponding term, \f$\log p_{outside}\f$ must not be added to the term.  For this reason, these ForceGenerators need to
/// be told if they are going to be used in a multiregion setting.  This is specified upon the construction of the 
/// ForceGenerators by a boolean argument to their constructor.
///
/// The algorithm in LSSEG that carries out the evolution of several level-set functions \f$\Phi_i\f$ in a multiresolution
/// setting is called \ref develop_multiregion_2D() "develop_multiregion_XD()".
///
/// \note Multiregion segmentation is a very unstable process and its very easy to get stuck in a local minimum when looking
/// for an optimal segmentation.  In general, this segmentation is quite sensitive to the initialization.  In order to obtain
/// a good multiregion segmentation of an image with LSSEG, the user should be prepared to experiment around a bit with various
/// settings in order to get an acceptable result.
///
///
/// \section section_Parts How the parts fit together
///
/// Here, we will attempt to describe an overview of the segmentation algorithm found in the
/// \ref sample_programs_page "sample program" called \ref segmentation2D.C "segmentation2D".
/// 
/// \subsection subsection_Initialization initialization of level-set functions
/// First of all, we define some initial partition of the image, from which to start the region development.  This is
/// carried out by a local function, \c initialize_levelsets(), which selects the use of one of the LSSEG functions:
/// \li \ref sphere() (two-region segmentation - initializes the level-set function so that its \em inside
/// region is a circle centered on the middle of the image)
/// \li \ref horizontal_sinusoidal_bands() (two-region segmentation - initializes the level-set function so that its inside
/// region becomes a set of horizontal bands across the image.
/// \li \ref multiregion_bands() (multiregion segmentation - initializes the set of level-set functions so that their \em
/// inside regions together becomes a set of horizontal bands across the image.
/// \li \ref random_scattered_voronoi() (multiregion segmentation - initializes the set of level-set functions so that their
/// \em inside regions together becomes a random, scattered partition of irregular fragments covering the image domain.
///
/// \subsection subsection_Multiresolution multiresolution
/// The first thing to mention is that we do not start a segmentation on the complete image right away.  Instead, we downscale
/// the image several times, run the segmentation on the smallest version on the image first, and then successively run
/// the segmentation on larger versions using the previous result as an initialization.  The reason for this approach is that
/// our PDE-based segmentation approach is computationally very expensive and do not converge very rapidly on fine grids
/// (i.e., on high-resolution images).  A segmentation process with several levels off resolution allows us to converge rapidly
///  on a coarse image, and hopefully get an initialization that is much closer to the desired result when running on higher
/// resolutions.
/// \image html multilevel.png "an illustration of the multilevel segmentation approach"
///
/// \subsection subsection_UseOfMask use of an image mask
/// It is possible for the user to specify a \ref Mask, which will limit the domain of the image on which the segmentation
/// will be carried out.  This can be useful in cases where only a specific part of the image is of interest for segmentation, 
/// and that the rest of the image will only "confuse" the segmentation process.  A mask is an single-channeled image 
/// whose pixels are either zero (off) or nonzero (on).  As such, it can be stored in any image format, as long as it is 
/// limited to a single channel.
///
/// \subsection subsection_SettingUpRegions "setting up regions"
/// 
/// The information needed about each region is:
///
/// \li the level-set function that describes its region
/// \li the scalar parameter describing the smoothness of its boundary
/// \li the \ref lsseg::ForceGenerator "ForceGenerator" that determines the
///  \ref anchor_NormalDirectionFlow "normal direction flow" components of its evolution.
///
/// These three pieces of information are bundled together in \ref lsseg::Region "Region" objects.
///
/// \subsection subsection_AssemblingMultichannelImage assembling the information
/// 
/// On each level in the multiresolution process, the program assembles a "working image", whose channels contain all the 
/// information that should steer the segmentation process on that level.  This may or may not include color channels, 
/// intensity channel or synthetized images ( \ref anchor_StructureTensor "structure tensor", 
/// \ref anchor_ScaleMeasure "scale measure") generated from the original image in order to enhance region differences. 
/// In the algorithm, this working image is called \c multichan in \ref segmentation2D.C "the code".  The channels going
///  into the assembly of this image can be specified by the user on the command line.
///
/// \subsection subsection_Prefiltering prefiltering
/// Before starting the segmentation process on the assembled image \c multichan, it is smoothed using an 
/// \ref anchor_ASmooth "anisotropic smoothing process" that preserves edges in the image.  The benefit of this smoothing 
/// is to remove fluctuations in order for the segmentation process to converge faster and prevent getting stuck in local
///  minima.  It also distributes more evenly the information found in the \ref anchor_StructureTensor "structure tensor" and 
/// \ref anchor_ScaleMeasure "scale measure" channels.  The level of smoothing to be used can be set from the command line.  
/// The user specifies a number, and when the
/// average pixel variation between two smoothing iterations is less that that number, the smoothing process stops.  For 
/// this reason, <em> smaller numbers lead to more smoothing </em>.  In general, smoothing is good, but too much smoothing
/// may remove important cues for the segmentation.  In addition, the smoothing is computationally expensive.
///
/// \subsection subsection_DevelopRegionXD call to the actual segmentation process
/// When the \ref lsseg::Region "regions" have been set up, and the working image (\c multichan) for a given level of resolution
/// has been assembled and smoothed, the segmentation process is started.  This is done through calls to the 
/// \ref develop_single_region_2D() function or the \ref develop_multiregion_2D() function, depending on whether we segment into
/// two or more regions.
///
/// \subsection subsection_PeriodicReinitialization periodic reinitialization
/// With regular intervals, the segmentation process is interrupted and the level-set functions
/// \ref anchor_Reinitialization "reinitialized" to the \ref anchor_SignedDistanceFunction "signed distance function" for the
/// region they currently describe.  This prevents situations where progress in the evolution of regions is  "blocked" due
/// to degenerate/extreme values in the gradient of the level-set functions.  (When reinitialized to the distance function,
/// the gradient becomes unity almost everywhere.
///
/// \subsection subsection_VisualizationOfResults visualization of results
/// Each time the segmentation process is interrupted for \ref anchor_Reinitialization "reinitialization" purposes, the algorithm
/// also updates a visualization window showing the progress of the segmentation (unless told not to do so).  
/// This is done by a call to the function \ref visualize_level_set() or \ref visualize_multisets(), which produces images
/// that fills the segmented region(s) with a color, and enhance the contour of the region.  The produced image is then
///  \ref lsseg::Image::superpose() "superposed" on the original image and shown in the visualization window.
///
/// Below you find the segmentation result of an image of a guppy.  The first segmentation contain two regions, the second
/// contain four.  The parameters for the four-region image had to be tuned more carefully.  For a higher resolution image,
/// click <a href="../images/segmented_guppy.png" target="_blank">here</a>.
/// \image html segmented_guppy_small.png "Left: original image. Center: segmented into two regions.Right: segmented into four regions"
/// The command line for the two-regions segmentation was <c>segmentation2D guppy.png -init 10</c>\n
/// The command line for the four-region segmentation was <c>segmentation2D guppy.png -init 4 -nr 4 -s 0.1 -lr 60 -i1 10000 </c>\n 
/// The multi-region algorithm had to be tuned much more carefully.
///
//=================== PAGE ON HOW LSSEG SMOOTHING WORKS ==============================
///
/// \page page_HowSmoothingWorks How Anisotropic Smoothing in LSSEG Works
///
/// \note Most of the information in this chapter is based on
/// \ref anchor_Tschumperle02 "the doctoral thesis of D. Tschumperle" and his 2006 paper 
/// \ref anchor_Tschumperle06 "Fast Anisotropic Smoothing of Multi-Valued Images using Curvature-Preserving PDE's".
/// This page is no substitute for reading these papers.  Rather, it is intended to give an overview that will
/// prepare the reader on the topic so that it might become easier to absorb the material presented in these papers.
/// 
/// \section section_SmoothingPDE Using partial differential equations for smoothing digital images
/// The reasons why we would want to smooth a digital image might be to remove noise, simplify the image or
/// remove local variations that might become obstacles for further processing (like segmentation).  
/// When using a partial differential equation (PDE) to smooth an image, we consider the image a discretized 
/// function over its domain, and apply a time-dependent PDE that makes the image evolve timestep by timestep
/// towards progressively simplified and smooth versions of itself.  The limit of such a process is generally
/// a image with a constant value everywhere; not a very interesting solution, but we usually stop the process
/// long before this happens.
/// All PDE's that we discuss here are on the form:
/// \f[\frac{\partial I}{\partial t} = R\f]
/// Here, \f$I\f$ represents the image, while \f$R\f$ represents an (image-dependent) smoothing term to be
/// further specified.
/// \section section_ObtainingEquation How do we obtain the equations
/// The simplest smoothing schemes can be derived directly from a \ref anchor_Functional "functional" using
/// the \ref anchor_EulerLagrange "Euler-Lagrange equation".  The functional in question gives a value that 
/// in some sense represents the amount of variation in the image, and the process of minimizing this functional
/// leads to a gradual transformation of the original image into simplified versions of itself.  We give here
/// two examples:
/// \subsection subsection_Tikhonov Example 1: isotropic smoothing
/// This is perhaps the simplest case one can think of, and it leads to a smoothing process that is specified
/// by the \ref anchor_HeatEquation "heat equation", and leads to linear diffusion.
/// The functional that we aim to minimize is the squared image gradient norm integrated over the image
/// domain \f$\Omega\f$:
/// \f[E(I) = \int_\Omega ||\nabla I||^2 d\Omega\f]
/// Which, through the use of the \ref anchor_EulerLagrange "Euler-Lagrange equation" gives us the following
/// PDE (the \ref anchor_HeatEquation "heat equation"):
/// \f{eqnarray*}  I_{(t=0)} & = & I_0 \\ \frac{\partial I}{\partial t} & = & \Delta I \f}
/// This smoothing process is \em isotropic; it smooths with equal strength in all spatial directions.  This
/// means that edges in the image will soon become blurred.  In LSSEG, we can apply this smoothing
/// scheme by calling the \ref nonlinear_gauss_filter_2D "nonlinear_gauss_filter_XD()" function and setting
/// the parameter \c p equal to zero.
/// \image html tikhonov.png "Example of progressive smoothing of an image using the heat equation"
/// \subsection subsection_TV Example 2: total variation flow
/// This time, we depart from the following functional:
/// \f[E(I) = \int_\Omega ||\nabla I|| d\Omega\f]
/// We see that it is somewhat similar to the integral in our last example, except that we have not squared 
/// the norm.  This functional measures the <em> total variation </em> of the image, and the PDE we obtain
/// from the \ref anchor_EulerLagrange "Euler-Lagrange equation" describes a process we call <em> total variation
/// flow</em>.  This process has several interesting \em theoretical properties; readers are refered to 
/// section 2.1.1 of \ref anchor_Brox05 "Brox' thesis" for an overview of these, and to section 2.3.2 to 2.3.4
/// of the same thesis for an examination of the resulting PDE.  As for the previous example, the limit of the
/// smoothing process is a constant image (and we will arrive there in finite time).  However, the intermediate 
/// images are quite different.  The edges in the image are well preserved for a long time, regions gradually 
/// merge with each other, and the intermediate images take on a segmentation-like appearance.
/// The PDE system is:
/// \f{eqnarray*}  I_{(t=0)} & = & I_0 \\ \frac{\partial I}{\partial t} & = & div\left(\frac{\nabla I}{||\nabla I||}\right)\f}
/// In LSSEG, we can apply total variation flow on an image using the function 
/// \ref nonlinear_gauss_filter_2D "nonlinear_gauss_filter_XD()" with the parameter \c p equal to 1.
/// \image html tvflow.png "Example of progressive smoothing of an image using total variation flow"
///
/// However, smoothing PDEs are not limited to those obtained by the
///  \ref anchor_EulerLagrange "Euler-Lagrange equations" on some functional.  Although this kind of model has
/// the teoretical benefit of reducing some explicitely defined criterion, it does not always have the effect
/// on the image that we would like.  In general, what we seek are schemes that, in a predictable way, smooth
/// the image regions with different strength depending on the amount of local variation, and with different
/// strength along the local image \em gradient and \em isophote directions.  For this, we can make modified
/// PDEs that no longer represent the gradient of some functional, and whose evolution does not lead to 
/// a minimization of some defined quantity, but which may in turn produce results that are tailored to our 
/// requirements.
/// \section section_Anisotropic Anisotropic smoothing formulations
/// If we have a look at the time-dependent evolution equation from the two examples above, we see that they
/// both can be written on the general form:
/// \f[\frac{\partial I}{\partial t} = div\left(g(||\nabla I||) \nabla I\right) \f]
/// with \f$g : R \rightarrow R\f$ a scalar function that is 1 in the first example and \f$\frac{1}{||\nabla I||}\f$
/// in the second.  Without departing from a functional, we can create customized schemes by specifying this
/// function \f$g\f$ directly.  What we want is to slow down diffusion in regions with high gradients (i.e., edges)
/// while still allowing diffusion to take place at normal speed in "flat" regions.
/// This idea was put forward by Perona and Malik in 1990, where they proposed the following choice of \f$g\f$:
/// \f[g(||\nabla I||) = \exp\left(-\frac{||\nabla I||^2}{K^2}\right)\f]
/// \f$K\f$ is here a fixed threshold that differentiates homogeneous regions from regions of contours.
/// Tschumperle gives a good theoretical explanation of this in section 2.1.2 of
/// \ref anchor_Tschumperle02 "his thesis", where he also shows that not only does this process lead to smoothing
/// whose strength depend on the local gradient; it also smoothes with different strength (which can be explicitely
/// determined) along the local image gradient and isophote.
/// \subsection subsection_DivergenceBasedPDEs Divergence-based PDEs
/// We can go further, however.  Above, \f$g\f$ is a scalar function that only depends on the image gradient norm.
/// A more general formulation would be:
/// \f[\frac{\partial I}{\partial t} = div\left(D \nabla I\right)\f]
/// where \f$ D: \Omega \rightarrow P(2)\f$ is a tensor field defined on the domain \f$\Omega\f$ and with values
/// in the space \f$P(2)\f$of tensors described by positive-definite 2x2-matrices.  (In case of three-dimensional
/// images, i.e., volumes, we use \f$P(3)\f$ instead of \f$P(2)\f$).  We call \f$D\f$ a field of 
/// \ref anchor_DiffusionTensor "diffusion tensors".  Contrary to \f$g\f$ described above, \f$D\f$ depends not on the
/// <em>image gradient</em> itself, but rather on the <em> position in the image </em>.  Moreover, as \f$D(x,y)\f$
/// is a tensor rather than a scalar, its effect on the smoothing term also varies with the local direction
/// of \f$\nabla I\f$ (matrix-vector multiplication \f$D \nabla I\f$).  The smoothing at position \f$(x, y)\f$
/// will in general be strongest along the direction specified by the eigenvector with largest eigenvalue of
/// \f$D(x,y)\f$, and weakest in the direction specified by the other eigenvector (which, due to the symmetry of
/// the matrix, is perpendicular on the first).  See the section on 
/// \ref section_SmoothingGeometry "definition of a smoothing geometry" for more on how the tensor field 
/// \f$D\f$ can be specified.  For more on divergence-based PDEs, refer to section 2.1.4 of 
/// \ref anchor_Tschumperle02 "Tschumperle's thesis".
/// \subsection subsection_Oriented1DLaplacians Smoothing based on oriented 1D Laplacians 
/// When expanding the expression on the right hand side of the equation above, the placement of \f$D\f$ inside
/// the \f$div\f$ term leads to a term which depend on the derivative of \f$D\f$.  For this reason, we cannot
/// in general say that the smoothing process locally follows the directions specified by the tensor field \f$D\f$.
/// D. Tschumperle shows that a more "predictable" process can be obtained by using the smoothing geometry
/// in a trace-based formulation:
/// \f[\frac{\partial I}{\partial t} = trace(TH) \f]
/// Here, \f$T\f$ is the smoothing geometry and \f$H = H(x, y)\f$ is the local \ref anchor_Hessian "Hessian matrix"
/// of the image.  For more details on this equation, and a discussion about how it differs from the divergence-
/// based formulation, please refer to  2.1.5 and section 2.4.1 respectively of 
/// \ref anchor_Tschumperle02 "Tschumperle's PhD-thesis".  In section 3.2.1 of the same thesis, Tschumperle 
/// shows that this equation has an analytical solution which is:
/// \f[I(t) = I_0 * G^{(T, t)}\f]
/// Where \f$I_0\f$ is the original image, \f$I(t)\f$ is the image at time \f$t\f$ of the evolution process, 
/// and \f$ G^{(T, t)}\f$ is a bivariate oriented Gaussian kernel:
/// \f[ G^{(T, t)}(x) = \frac{1}{4\pi t}\exp\left(-\frac{[x,y]T^{-1}[x,y]^T}{4t}\right)\f]
///
/// \section section_SmoothingGeometry Defining a smoothing geometry
/// The purpose of the tensor field referred to as the "smoothing geometry" is to define local smothing directions
/// at each point \f$(x, y)\f$ in the image.  Typically, in regular areas of the image, we want to smooth 
/// equally in all directions, whereas in regions with important gradients, it is important to smooth the image
/// in a way that does not destroy the gradient, i.e., smooth much less \em along the gradient than in the 
/// direction \em perpendicular to it (the isophote direction).  Thus, in order to define the smoothing geometry
/// tensor field, we should, for each point in the image, define two directions (2D-vectors) and the smoothing
/// strength along these who directions (scalars).  In order to define these, we base ourselves on a smoothed
/// version of the image's \ref anchor_StructureTensor "structure tensor".  (The reason we smooth the structure 
/// tensor is to determine a gradient for each pixel that is somewhat less "local", thus representing structures
/// in the image that go beyond the size of one pixel and also obtaining gradient estimates that are more robust
/// against image noise).
///
/// For a given image point \f$(x, y)\f$, the eigenvectors of the structure tensor represent the two privileged
/// smoothing directions.  Due to the smoothing of the structure tensor, we cannot directly refer to these
/// eigenvectors as the gradient and isophote of the image at \f$(x, y)\f$, but they are perpendicular and 
/// represent the direction along which the image varies the most (and thus were we should smooth little) and
/// the direction along which the image varies the least (and thus where we should smooth more).  The corresponding
/// eigenvalues give information about strong the variation is.  If the two eigenvalues are identical, then the
/// image variation around this point is more or less the same in all directions, and the smoothing geometry
/// should be \em isotropic here.  If the two eigenvalues are different, we should smooth more along the direction
/// with the smaller eigenvalue than along the direction with the bigger one, and the difference in smoothing strength
/// (the anisotropy) should be reflected by the difference between these two values.
/// If we name the two directions \f$u\f$ and \f$v\f$, and choose two corresponding smoothing strenths \f$c_1\f$
/// and \f$c_2\f$, the smoothing geometry matrix representation for a given image point becomes:
/// \f[T = c_1 u u^T + c_2 v v^T\f]
/// In LSSEG, we have a function called \ref compute_smoothing_geometry_2D() "compute_smoothing_geometry_XD()"
/// that computes a smoothing geometry based on an image's structure tensor.  Here, the two smoothing directions
/// \f$u\f$ and \f$v\f$ are the eigenvectors of the structure tensor, and the corresponding smoothing strengths
/// are taken to be: 
/// \f[c_1 = \frac{1}{(1 + \lambda_u + \lambda_v)^{p_1}}\;,\;\;\;\;c_2 = \frac{1}{(1 + \lambda_u + \lambda_v)^{p_2}} \f]
/// where \f$\lambda_u\f$ and \f$\lambda_v\f$ are the two eigenvalues of the structure tensor corresponding to 
/// \f$u\f$ and \f$v\f$ respectively.  This choice is discussed briefly in section 2.1 of 
/// \ref anchor_Tschumperle06 "[Tschumperle06]".  Note that as long as \f$p_1 < p_2\f$, then the smoothing will 
/// be stronger along direction \f$u\f$ than along direction \f$v\f$, and this difference will be more pronounced
/// as the biggest eigenvalue (an thus the denominator of the fractions above) grows.
/// \note Everything said above is directly generalizable to 3D-images (volumes), and LSSEG handles 2D and 3D 
/// smoothing in the same way.
///
/// \section section_MultichannelImages Images with multiple channels
/// When smoothing images containing multiple \ref anchor_Channel "channels", it is usually beneficial to smooth
/// all image channels <em>using the same smoothing geometry</em>.  If we smooth each image channel separately
/// using a smoothing geometry computed independently of the other image channels, we risk ending up with an image
/// where image edges are positioned differently in each channel, giving a undersirable blurred effect.
/// Instead, we compute a common smoothing geometry based on the multi-channel definition of the image's
/// \ref anchor_StructureTensor "structure tensor".  This is covered in detail in section 2.2 of 
/// \ref anchor_Tschumperle02 "Tschumperle's PhD thesis".
///
/// \section section_CurvaturePreserving Smoothing that preserves curvatures
/// We saw earlier that the trace-based formulation behaves locally as an oriented gaussian smoothing.  The
/// strength and orientation of this smoothing is directly reflected by the given structure tensor \f$T\f$.
/// This does a good job in preserving edges, but on \em curved structures, especially sharp corners, this
/// leads to an undesirable rounding effect.  The problem is that the gaussian kernel does not take curvatures
/// into account.  As a remedy, Tschumperle proposes an approach which can be seen as local filtering of the
/// image with <em> curved gaussian kernels </em>.  He develops the theory and algorithm in section 3 and 4 of 
/// his paper 
/// \ref anchor_Tschumperle06 "Fast Anisotropic Smoothing of Multi-Valued Images using Curvature-Preserving PDE's".
/// In short, he splits the smoothing geometry <em> tensor field</em> \f$T\f$ into a sum of <em>vector fields</em>, 
/// and then performs \ref anchor_LIC "line integral convolutions" of the image with each of these vector
/// fields and a \em univariate gaussian kernel function: 
/// \f[G_t(s) = \frac{1}{4t}\exp\left(-\frac{s^2}{4t}\right)\f]
/// \f$t\f$ is here the timestep used.  In other words he carries out 1-dimensional convolutions of the image
/// with \f$G_t(s)\f$ <em> along the streamlines defined by the vector fields obtained from \f$T\f$</em>.
/// \subsection subsection_TensorSplitting decomposing the tensor field
/// In order to define vector fields along which to carry out \ref anchor_LIC "line integral convolution",
/// Tschumperle decompose the smoothing geometry \f$T\f$ in the following way:
/// \f[T = \frac{2}{\pi} \int_{\alpha = 0}^\pi (\sqrt{T}a_\alpha)(\sqrt{T}a_\alpha)^T\f]
/// where \f$a_\alpha = [\cos \alpha, \sin \alpha]^T\f$ and \f$\sqrt{T}\f$ is a tensor field that is the
/// "square root" of \f$T\f$, in the sense that \f$(\sqrt{T})^2=T\f$.  If \f$T = c_1 uu^T + c_2 vv^T\f$, then
///  \f$\sqrt{T} = \sqrt{c_1} uu^T + \sqrt{c_2} vv^T\f$.
/// \subsection subsection_ResultingPDE The resulting partial differential equation
/// Formally, the PDE that Tschumperle solves in the scheme proposed is:
/// \f[\frac{\partial I}{\partial t} = trace(TH) + \frac{2}{\pi}\nabla I^T \int_{\alpha=0}^\pi J_{\sqrt{T}a_\alpha}\sqrt{T}a_\alpha d\alpha\f]
/// where \f$J_{\sqrt{T}a_\alpha}\f$ stands for the Jacobian of the vector field
/// \f$\Omega \rightarrow \sqrt{T}a_\alpha\f$.
/// However, the system is not solved using the "standard" method of discrete differencing on a regular grid.
/// Instead, he demonstrates how this can be interpreted as a combination of 1D heat flows along the streamlines of
/// the various vector fields \f$\sqrt{T}a_\alpha\f$.
/// \section section_GreycStoration The LSSEG implementation of Tschumperles curvature-preserving smoothing algorithm
/// We have implemented a 2D/3D version of the algorithm outlined above, which can be found in the 
/// example program \ref smoothcurvpres.C "app/smoothcurvpres.C" that comes with LSSEG.  Here, we will briefly 
/// present what is going on inside this implementation.
/// -# First, the image is read from file.  This can be a normal 2D image in any standard format, or alternatively
///    a \em stack representing a 3D image.  Stacks can be made with the program
///    \ref generate_stack.C "app/generate_stack.C"
/// -# Then we enter a loop on the number of timesteps used.  Note that since we are using
///    \ref anchor_LIC "line integral convolution" to compute the changes, each timestep can be arbitrarily 
///    large without affecting the stability (the nonlinearity \em will be affected, though).
/// -# Depending on whether the image in question is 2D or 3D, the 2D/3D
///    \ref anchor_StructureTensor "structure tensor" is now computed.  Note that the image has been blurred
///    a little prior to this.
/// -# Based on the structure tensor, the smoothing geometry is now computed, by a call to 
///    \ref compute_smoothing_geometry_2D "compute_smoothing_geometry_2D()" or 
///    \ref compute_smoothing_geometry_3D "compute_smoothing_geometry_3D()", depending on whether the image is
///    two or three dimensional.
/// -# Now we carry out the splitting of the smoothing geometry tensor field into a sum of vector fields, 
///    and for each such field we carry out the actual \ref anchor_LIC "line integral convolution".
/// -# Finally, the result is displayed and, in case of 3D images, the result is saved as a \em stack to file.
///
/// \note due to problems with the functions \ref LIC_2D() and \ref LIC_3D(), these are substituted with the
/// less precise functions \ref LIC_2D_FS() and \ref LIC_3D_FS().  Read the note in the documentation of these
/// function for an explanation.
/// \image html lic.png "image smoothed using our implementation of Tschumperle's algorithm.  The parameters are the default ones, and the number of timesteps are 1, 4 and 8 (original image to the left)"

