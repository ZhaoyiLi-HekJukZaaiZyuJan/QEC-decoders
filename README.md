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
  <h2>Testing the Cppflow Package</h2>
  <p>
  To test the Cppflow package, go to <code>/CppflowTest</code> and run the following code to compile <code>main.cpp</code>:
  <pre>
  g++ -std=c++17 -o main.out -I ../src/libs/cppflow-master/include/ main.cpp -ltensorflow
  </pre>
  Running the code ./main.out would generate the following output:
  <pre>
  (tensor: shape=[3], dtype=TF_DOUBLE, data=[2 3 4])
  </pre>
  <p>
</section>

<section>
<h2>Training a model</h2>

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
<p>This program accepts command-line arguments, which can be used to customize its behavior. The available options for each subprogram are listed in the following table, along with any options that are not available. </p>
<table>
  <tr>
    <th>Option</th>
    <th>Description</th>
    <th>Default Value</th>
    <th>ML3D</th>
    <th>ML2D</th>
  </tr>
  <tr>
    <td>-f, --fname</td>
    <td>Output filename</td>
    <td>N/A</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
  <tr>
    <td>-d, --directory</td>
    <td>Model directory, enter the directory that contains the folder /models</td>
    <td>../src</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
  <tr>
    <td>-m, --model</td>
    <td>Model name</td>
    <td>model,L=5(7),layer=5x512,epochs=1000,p=</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
  <tr>
    <td>-s, --surf_type</td>
    <td>Surface type</td>
    <td>subPLANE</td>
    <td>&#x2705</td>
    <td>subTORUS</td>
  </tr>
  <tr>
    <td>--lmin</td>
    <td>Minimal size of lattice</td>
    <td>3</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
  <tr>
    <td>--lmax</td>
    <td>Maximal size of lattice</td>
    <td>17</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
  <tr>
    <td>-l</td>
    <td>Level</td>
    <td>0</td>
    <td>&#x2705</td>
    <td></td>
  </tr>
  <tr>
    <td>-n</td>
    <td>Number of trials</td>
    <td>10000</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
  <tr>
    <td>-v</td>
    <td>Verbosity switch which can take value 0-2</td>
    <td>0</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
  <tr>
    <td>--generate</td>
    <td>Generation mode switch. Do not use this option unless you want to generate training data</td>
    <td>FALSE</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
  <tr>
    <td>--Np</td>
    <td>Z error p points</td>
    <td>10</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
  <tr>
    <td>--pmin</td>
    <td>Minimal Z error probability</td>
    <td>0.001</td>
    <td>&#x2705</td>
    <td>0.01</td>
  </tr>
  <tr>
    <td>--pmax</td>
    <td>MaximalZ error probability</td>
<td>0.008</td>
<td>&#x2705</td>
<td>0</td>
  </tr>
  <tr>
    <td>-c, --code_model</td>
    <td>Code model</td>
    <td>2D</td>
    <td>&#x2705</td>
    <td></td>
  </tr>
  <tr>
    <td>--test</td>
    <td>Test switch</td>
    <td>0</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
  <tr>
    <td>--sweep</td>
    <td>Sweep switch</td>
    <td>0</td>
    <td>&#x2705</td>
    <td></td>
  </tr>
  <tr>
    <td>-N, --noise_model</td>
    <td>Noise model</td>
    <td>INDEP</td>
    <td>&#x2705</td>
    <td>DEPOL</td>
  </tr>
  <tr>
    <td>--seed</td>
    <td>Seed switch</td>
    <td>0</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
  <tr>
    <td>--thread</td>
    <td>Thread switch</td>
    <td>0</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
  <tr>
    <td>--use_env</td>
    <td>Use environment variables</td>
    <td>0</td>
    <td>&#x2705</td>
    <td>&#x2705</td>
  </tr>
</table>

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


<h2>Examples</h2>

Here are some examples of how to use the program with different command line arguments:

<pre>
./simulate -s subTORUS --pmin 0 --pmax 0.18  --Np 25 -n 1000 --Lmin 5 --Lmax 5 -v 1 -d ~/ML/ -m "model,L=5(7),layer=3x128,epochs=10000,p=" --decode_with_NN
</pre>
rcx
This command will run simulations for a subTORUS surface with a range of error probabilities between 0 and 0.18, which passes the transition error rate of around 0.15, using 25 p points and 1000 trials per point. The size of the mesh will be fixed at 5, and a neural network model with architecture 3x128 and 10000 epochs will be used for decoding. The results will be saved to the specified output directory with a filename that includes the model and input parameters. [expand a bit on the directory ]

<pre>
./simulate -s subTORUS --pmin 0 --pmax 0.12  --Np 25 -n 1000 --Lmin 3 --Lmax 20 -v 1 -d ~/ML/  --fname test.out
</pre>
<pre>
./simulate -s subTORUS --pmin 0.036 --Np 10 --Lmin 10 -v 1 --test --make_corrections -d /scratch/users/ladmon/ML/ -m "model_h,L=5(7),layer=3x128,epochs=100000,p=0.036" --binary
</pre>
<pre>
./simulate -s subTORUS --pmin 0.036 --Np 10 --Lmin 10 -v 1 --test --make_corrections -m "model,L=5(7),layer=5x512,epochs=1000,p=0.04" --binary
</pre>
<pre>
./simulate -s subTORUS --pmin 0.02 --pmax 0.02  --Np 20 -n 1 --Lmin 7 -v 1 --generate -d ~/ML
</pre>