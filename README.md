<!--
*** Thanks for checking out this README Template. If you have a suggestion that would
*** make this better, please fork the repo and create a pull request or simply open
*** an issue with the tag "enhancement".
*** Thanks again! Now go create something AMAZING! :D
-->





<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
![code-lines-shield]



<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/StudDavid/biofilm_growth_modeling">
  </a>

  <h3 align="center">iGEM 2020 project: Biofilm growth simulation</h3>

  <p align="center">
    This is repository contains the code for simulating the growth of a Biofilm after attachment using Molecular dynamics simulation methods. The project was build as part of the iGEM 2020 Contest. The numerical model and simulation were implmented by the Team Hannover. Functionalities to make the software tool more accesible were added in collaboratio with the Team TU Darmstadt. The BiofilmSimulation software can be accessd a a python package, documented under https://pypi.org/project/BiofilmSimulation/1.0.0/.
    <br />
    <a href="https://github.com/StudDavid/biofilm_growth_modeling"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/StudDavid/biofilm_growth_modeling">View Demo</a>
    ·
    <a href="https://github.com/StudDavid/biofilm_growth_modeling/issues">Report Bug</a>
    ·
    <a href=https://github.com/StudDavid/biofilm_growth_modeling/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
## Table of Contents

* [About the Project](#about-the-project)
  * [Built With](#built-with)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
* [Usage](#usage)
* [Roadmap](#roadmap)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)
* [iGEM Competition](#igem-competition)
* [Acknowledgements](#acknowledgements)



<!-- ABOUT THE PROJECT -->
## About The Project

The project aims to simulate the growth of a biofilm in early stages. An biofilm is an consortium of bacteria embedded in a extracelluar matrix consisting of EPS  (extracellular polymeric substances). Origin of such biofilms is the attachment of initial bacteria to a surface.

As part of the iGEM Competition 2020, the Hannover Team designed a sensor based on biological cell, which is capable of detecting the adhernce of a biofilm to at an early stage. The sensor can be attached to implant surfaces and used as a diagnostical tool.
Therefore, we are interested in the growth behaviour of biofilms in an early stage.
We use computational methods of Molecular Dynamics simulation and a biophysical approach to model the biofilm growth. Furthermore the BiofilmSimulation module provides functions for visualization of the data generated in the simulation.

### Built With

* [SciPy](https://www.scipy.org/)
* [NumPy](https://numpy.org/)
* [Pandas](https://pandas.pydata.org/)


<!-- GETTING STARTED -->
## Getting Started

This is an example of how you may set up your project locally.
To get a local copy up and running follow these simple example steps.

### Prerequisites 
We recommend installing anaconda on your machine. Anaconda provides many functionalities including an easy way to set up python enviroments. 
Check out the [Anaconda installation instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Installation

You will need a few python packages to run the simulation on your local machine. You can eithercreate a conda enviroment from the 
`conda_env.yml` file or download the packages on your own via pip.
We provide step-by-step instructions on how to get our software running using a conda enviroment.

1. Clone the repo
```sh
git clone https://github.com/StudDavid/biofilm_growth_modeling.git
```
2. Navigate in the root folder (the folder which contains the `setup.py`) and run

```sh
conda env create -f env/iGEM-biofilm-model.yml
```

3. Start the enviroment by running 
```sh
conda activate iGEM-biofilm-model
```

4. Install our software tool as a python package in your conda enviroment with pip (package installer for Python). The package is listed at https://pypi.org/project/BiofilmSimulation/1.0.0/ and can be installed by running 

```sh
pip install BiofilmSimulation
```

Another option is to run
```sh
pip install .
```
in the root folder. This will run the setup.py file and also adds our BiofilmSimulation module to your enviroment.
In both cases, if the installation was succesfull the anaconda prompt will print 
```sh
Successfully built BiofilmSimulation
Installing collected packages: BiofilmSimulation
Successfully installed BiofilmSimulation-1.0.0
```
To be able to save the 2D animation of the simulation, you also need the `ffmpeg - Writer`. Install it in your enviromnet with:
```sh
conda install -c conda-forge ffmpeg
```

<!-- USAGE EXAMPLES -->
## Usage
You can also directly start a simulation by running the `example_usage.py` script. 
Activate the `iGEM-biofilm-model` enviroment as described above and navigate to the root folder. 
Type 
```sh
python example/example_usage.py
```
This will start a dialog in the console in which you can specify the bacteria strain, the number of initial bacteria and the simulation duration.

Another usage example is provided in form of a jupyter notebook. Make sure to check it out to see the functionalities provided. 
If you are using anaconda, jupyter will already be installed. To start the example jupyter notebook run

```sh
jupyter-notebook
```
Then connect to 
```sh
http://localhost:8888/
```
and navigate to the folder, in which you cloned the repository. Start the notebook by clicking on the `example.ipynb` file.



<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/StudDavid/biofilm_growth_modeling/issues) for a list of proposed features (and known issues).


<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

* David Theidel - theidel AT stud dot uni-hannover dot de

* Project Link: [https://github.com/StudDavid/biofilm_growth_modeling](https://github.com/StudDavid/biofilm_growth_modeling)

## iGEM Competition
* To find out more abour the awesome **iGEM Comeptition** check our [https://2020.igem.org/](https://2020.igem.org/)

* More about the project of the Hannover Team can be found here: [https://2020.igem.org/Team:Hannover](https://2020.igem.org/Team:Hannover)

* Details on the used methods: [https://2020.igem.org/Team:Hannover/Software](https://2020.igem.org/Team:Hannover/Software)

<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements
* [GitHub Emoji Cheat Sheet](https://www.webpagefx.com/tools/emoji-cheat-sheet)
* [Img Shields](https://shields.io)
* [Choose an Open Source License](https://choosealicense.com)
* [GitHub Pages](https://pages.github.com)
* [Animate.css](https://daneden.github.io/animate.css)
* [Loaders.css](https://connoratherton.com/loaders)
* [Slick Carousel](https://kenwheeler.github.io/slick)
* [Smooth Scroll](https://github.com/cferdinandi/smooth-scroll)
* [Sticky Kit](http://leafo.net/sticky-kit)
* [JVectorMap](http://jvectormap.com)
* [Font Awesome](https://fontawesome.com)





<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/othneildrew/Best-README-Template.svg?style=flat-square
[contributors-url]: https://github.com/StudDavid/biofilm_growth_modeling/contributors
[stars-shield]: https://img.shields.io/github/stars/StudDavid/biofilm_growth_modeling?style=social
[stars-url]: https://github.com/StudDavid/biofilm_growth_modeling/stargazers
[forks-shield]: https://img.shields.io/github/forks/othneildrew/Best-README-Template.svg?style=flat-square
[issues-shield]: https://img.shields.io/github/issues/othneildrew/Best-README-Template.svg?style=flat-square
[issues-url]: https://github.com/StudDavid/biofilm_growth_modeling/issues
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=flat-square
[license-url]: https://github.com/StudDavid/biofilm_growth_modeling/blob/master/LICENSE
[code-lines-shield]: https://img.shields.io/tokei/lines/github/StudDavid/biofilm_growth_modeling
