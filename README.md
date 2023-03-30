<h1>QEC-decoders</h1>
<p>This repository contains several decoders for topological quantum error correction. The decoders are implemented in C++ and use the Blossom V algorithm for maximum-weight perfect matching.</p>
<h2>Usage</h2>
<ol>
  <li>Clean the object files by running <code>make clean</code> in the <code>src/blossom5-v2.05.src</code> directory:</li>
    <pre>cd src/blossom5-v2.05.src
make clean</pre>
  <li>Clean the project files:</li>
    <pre>make clean</pre>
  <li>Make the project files:</li>
    <pre>make</pre>
  <li>Run the executable:</li>
    <pre>./simulate</pre>
  <li>Parameters for the simulations are passed in according to the comments inside the code.</li>
</ol>
<p>If you encounter any issues or have questions, please feel free to open an issue or submit a pull request.</p>


