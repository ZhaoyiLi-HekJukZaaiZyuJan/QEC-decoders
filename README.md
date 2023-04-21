<h1>QEC-decoders</h1>
<p>This repository contains several decoders for topological quantum error correction. The decoders are implemented in C++ and use the Blossom V algorithm for maximum-weight perfect matching.</p>
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
<p>If your system does not have TensorFlow installed, you can build TensorFlow with the following commands:</p>
<pre>FILENAME=<correct file name>
wget -q --no-check-certificate https://storage.googleapis.com/tensorflow/libtensorflow/${FILENAME}
sudo tar -C /usr/local -xzf ${FILENAME}</pre>
<p>Check <a href="https://www.tensorflow.org/install/lang_c">https://www.tensorflow.org/install/lang_c</a> for appropriate versions.</p>
<p>Try to compile the source in TF test to test TensorFlow.</p>
