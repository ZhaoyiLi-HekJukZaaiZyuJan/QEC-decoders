<h1>QEC-decoders</h1>
<p>This repository contains several decoders for topological quantum error correction. The decoders are implemented in C++ and use the Blossom V algorithm for maximum-weight perfect matching.</p>

<section>
  <h2>Setting up TensorFlow Dependencies</h2>
  <p>If your system does not have TensorFlow installed, you can build TensorFlow with the following commands:</p>
  <pre>
FILENAME= # Set your filename here, e.g. libtensorflow-cpu-linux-x86_64-2.11.0.tar.gz
wget -q --no-check-certificate https://storage.googleapis.com/tensorflow/libtensorflow/${FILENAME}
sudo tar -C /usr/local -xzf ${FILENAME}
  </pre>
  <p>Note: Be sure to check <a href="https://www.tensorflow.org/install/lang_c">https://www.tensorflow.org/install/lang_c</a> for appropriate versions.</p>
  <p>To ensure that your local TensorFlow dependencies are correctly set up, you can try compiling the source code in the TensorFlow test folder:</p>
  <pre>
cd tensorflow
./configure
bazel test tensorflow/test/...
  </pre>
</section>

<section>
<h2>Usage</h2>
<ol>
  <li>Make the src files:</li>
    <pre>cd src/libs/blossom5-v2.05.src <br>make <br>cd ../../.. </pre>
  <li>Make the project files:</li>
    <pre>make</pre>
  <li>Run the executable:</li>
    <pre>./simulate</pre>
  <li>Parameters for the simulations are passed in according to the comments inside the code.</li>
</ol>
<p>If you encounter any issues or have questions, please feel free to open an issue or submit a pull request.</p>
</section>

<section>
  <h2>Variables</h2>
<p>The following command line arguments can be passed to the program:</p>
<ul>
  <li>-f, --fname: Output filename.</li>
  <li>-s, --surf_type: Surface type. Default is subPLANE.</li>
  <li>--lmin: Minimal size of mesh. Default is 3.</li>
  <li>--lmax: Maximal size of mesh. Default is 17.</li>
  <li>-l: Level. Default is 0.</li>
  <li>-n: Number of trials. Default is 10000.</li>
  <li>-v: Verbosity switch. Default is 0.</li>
  <li>--Np: Z error p points. Default is 10.</li>
  <li>--pmin: Minimal Z error probability. Default is 0.001.</li>
  <li>--pmax: Maximal Z error probability. Default is 0.008.</li>
  <li>-c, --code_model: Code model. Default is 2D.</li>
  <li>--test: Test switch. Default is 0.</li>
  <li>--sweep: Sweep switch. Default is 0.</li>
  <li>-N, --noise_model: Noise model. Default is INDEP.</li>
  <li>--seed: Seed switch. Default is 0.</li>
  <li>--thread: Thread switch. Default is 0.</li>
  <li>--use_env: Use environment variables. Default is 0.</li>
</ul>
<p>To run the program with default settings, simply run:</p>
<pre>./simulate</pre>

<h2>Functions</h2>
<p>The following functions are defined in the code:</p>
<ul>
  <li><code>testDecoding</code>: This function tests the decoding process for a given sub-cluster.</li>
  <li><code>SweepTestDecoding</code>: This function performs a sweep test on the decoding process for a given sub-cluster.</li>
  <li><code>loopDecoding</code>: This function performs the decoding process for a given range of parameters and saves the results to a file.</li>
</ul>

<h2>Inputs and Outputs</h2>
<p>The inputs and outputs of the program are defined by the command line arguments and the output file. The output file contains the following data:</p>
<ul>
  <li>L: The size of the mesh.</li>
  <li>p_error: The probability of a Z error.</li>
  <li>num_success: The number of successful decoding attempts.</li>
</ul>
<p>If no output filename is specified, the program will create a default output file name based on the input parameters.</p>
</section>
