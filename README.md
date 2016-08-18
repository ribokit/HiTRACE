# HiTRACE (High Throughput Robust Analysis of Capillary Electrophoresis)

**HiTRACE** is a collection of *MATLAB* scripts to automate the key tasks in large-scale nucleic acid CE analysis, including the profile alignment that has heretofore been a rate-limiting step in the highest throughput experiments. It has been intensively used for quantitating data for RNA and DNA based on the mutate-and-map methodology, chromatin footprinting, and other high-throughput structure mapping techniques.

An online user-friendly GUI is available at the [**HiTRACE Web**](http://hitrace.org/).

## Installation

**HiTRACE** requires *MATLAB* version `>= R2011a` and `<= R2014a`. Later versions of *MATLAB* are incompatible with codes that involves interactive interface handling. For *Mac OS X* users with version `>= 10.10`, you may need this [patch](https://www.mathworks.com/support/bugreports/1098655) to settle a *Java* exception.

To install **HiTRACE**, simply:

- Download the zip or tar file of the repository and unpack; or 
```bash
git clone https://github.com/hitrace/HiTRACE.git
```

- In *MATLAB*, go to "**Set Path**". Then "**Add with Subfolders**" of the target `path/to/HiTRACE/Scripts/`.

## Usage 

Documentation (*MATLAB* tutorial) is available at https://hitrace.github.io/HiTRACE/.

## License

Copyright &copy; of **HiTRACE** _Source Code_ is described in [LICENSE.md](https://github.com/hitrace/HiTRACE/blob/master/LICENSE.md).

## Reference

>Lee, S., Kim, H., Tian, S., Lee, T., Yoon, S., and Das, R. (**2015**)<br/>
>[Automated band annotation for RNA structure probing experiments with numerous capillary electrophoresis profiles](http://bioinformatics.oxfordjournals.org/content/31/17/2808.abstract)<br/>
>*Bioinformatics* **31 (17)**: 2808 - 2815.

>Kladwang, W., Mann, T.H., Becka, A., Tian, S., Kim, H., Yoon, S., and Das, R. (**2014**)<br/>
>[Standardization of RNA chemical mapping experiments](http://pubs.acs.org/doi/abs/10.1021/bi5003426)<br/>
>*Biochemistry* **53 (19)**: 3063 - 3065.

>Kim, H., Cordero, P., Das, R., and Yoon, S. (**2013**)<br/>
>[HiTRACE-Web: an online tool for robust analysis of high-throughput capillary electrophoresis](http://nar.oxfordjournals.org/content/41/W1/W492)<br/>
>*Nucleic Acid Research* **41 (W1)**: W492 - W498.

>Yoon, S., Kim, J., Hum, J., Kim, H., Park, S., Kladwang, W., and Das, R. (**2011**)<br/>
>[HiTRACE: High-throughput robust analysis for capillary electrophoresis](http://bioinformatics.oxfordjournals.org/content/27/13/1798)<br/>
>*Bioinformatics* **27 (13)**: 1798 - 1805.

<br/>
Developed by **Das lab** (_Leland Stanford Junior University_) and **Yoon lab** (_Seoul National University_) and colleagues.
<br/>
README by [**t47**](http://t47.io/), *April 2016*.
