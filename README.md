<h1>QEC-decoders</h1>
<p>This repository contains several decoders for topological quantum error correction. The decoders are implemented in C++ and use the Blossom V algorithm for maximum-weight perfect matching.</p>

<section>
  <h2>Setting up TensorFlow Dependencies</h2>
  <p>If your system does not have TensorFlow installed, you can build TensorFlow with the following commands:</p>
  <pre>
FILENAME=#(Set your filename here, e.g. libtensorflow-cpu-linux-x86_64-2.15.0.tar.gz)
wget -q --no-check-certificate https://storage.googleapis.com/tensorflow/versions/(add versions)/libtensorflow/${FILENAME}
sudo tar -C /usr/local -xzf ${FILENAME}
export DYLD_LIBRARY_PATH=/usr/local/lib:$DYLD_LIBRARY_PATH
./hello_tf</pre>

  <p>To ensure that your local TensorFlow dependencies are correctly set up, you can try compiling the source code in the TFtest folder:</p>
  <pre>
gcc hello_tf.c -ltensorflow -o hello_tf
./hello_tf</pre>
  <h2>Testing the Cppflow Package</h2>
  <p>
  To test the Cppflow package, go to <code>CppflowTest</code> and run the following code to compile <code>main.cpp</code>:
  <pre>
  g++ -std=c++17 -o main.out -I ../src/libs/cppflow-master/include/ main.cpp -ltensorflow
  </pre>
  Running the code ./main.out would generate the following output:
  <pre>
  (tensor: shape=[3], dtype=TF_DOUBLE, data=[2 3 4])
  </pre>
  <p>
  <p>Note: Be sure to check <a href="https://www.tensorflow.org/install/lang_c">https://www.tensorflow.org/install/lang_c</a> for more details.</p>
</section>

<section>
<h2>Usage</h2>
<ol>
  <li>Make the src files:</li>
    <pre>cd src/libs/blossom5-v2.05.src <br>make <br>cd ../../.. </pre>
  <li>Go to a project folder, e.g. 2D, and make the project files:</li>
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
    <p>This flag enables testing mode, which will run a single decoding simulation for the given size of the mesh and a specified error probability. The simulation process will be printed to the console. The default value of this flag is 0, which means testing mode is disabled.</p>
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

Here’s an HTML block you can paste straight into your README.md (GitHub renders inline HTML):

<h3>Coordinate Convention</h3>

<p>
This code labels qubits/checks with a 4-tuple <code>(x, y, z, l)</code> on a 3D periodic grid.
The spatial sizes are <code>(L, M, N)</code> and the sublattice index is <code>l</code>. For qubits, <code>l</code> ranges from 0 to 2.
</p>

<table>
  <thead>
    <tr><th>Symbol</th><th>Meaning</th></tr>
  </thead>
  <tbody>
    <tr><td><code>x, y, z</code></td><td>Spatial coordinates on a torus of size <code>(L, M, N)</code></td></tr>
    <tr><td><code>l</code></td><td>Internal layer / sublattice index</td></tr>
    <tr><td><code>L, M, N</code></td><td>Lattice extents; provided via a size holder <code>S</code> with <code>S.x=L, S.y=M, S.z=N</code></td></tr>
  </tbody>
</table>

<p>
A site has a unique linear ID computed by <code>coord::hash</code>:
</p>

<pre><code>id = l·L·M·N  +  z·L·M  +  y·L  +  x
</code></pre>

<p>
The inverse mapping (flattened ID → coordinates) is implemented by
<code>coord::coord(const int&amp; c, const coord&amp; S)</code>, which unpacks with modular arithmetic:
</p>

<pre><code>x =  c % (M·L·N) % (M·L) % L
y =  c % (M·L·N) % (M·L) / L
z =  c % (M·L·N) / (M·L)
l =  c / (M·L·N)
</code></pre>

<p>
Dual-lattice face qubits are addressed by functions <code>getFaceQubits(...)</code>.
They select neighboring hashed coordinates on the appropriate face (xy / zx / zy),
with parity rules on <code>x, y, z</code> and wrap-around via <code>divmod</code>.
The compact overload <code>getFaceQubits(S, k)</code> maps a face label
<code>k ∈ {0,3,4,6,7,8,9,10,11}</code> to the corresponding neighbor’s linear ID.
</p>

<p>
For 2D substructures, <code>subcoord</code> uses <code>(x, y, l)</code> with sizes
<code>(L, M)</code> and analogous hash/unhash:
</p>

<pre><code>sub-id = l·L·M  +  y·L  +  x
x = c % (M·L) % L
y = c % (M·L) / L
l = c / (M·L)
</code></pre>

<p>
Example usage:
</p>

<pre><code>// Sizes: L=4, M=3, N=2
coord S(4, 3, 2, 0);          // size holder (l not used)
coord q(2, 1, 0, 1);          // a site
int id = q.hash(S);           // flatten → unique integer
coord r(id, S);               // unflatten → (x,y,z,l) == (2,1,0,1)
</code></pre>