# Human-in-the-Loop Parameter Selection for GFP Stack Segmentation (HITL)

This package contains Python code implementing a human-in-the-loop (HITL) approach for parameter selection in GFP (Green Fluorescent Protein) stack segmentation. 
The process refines segmentation parameters through multiple phases, incorporating expert feedback at each stage.

## Overview

The code facilitates iterative segmentation parameter selection / optimization. We use 4 phases. The number of phases can be chosen as required:

The first phase parameters were selected using estimated parameters by measuring structure sizes in imageJ<ref>.
For all phases, the following steps are followed:
1.  **Segmentation and Expert Review:** Image stacks are segmented using these parameters. Experts review the results and provide feedback on the quality of the segmentation.
2.  **Parameter Refinement:** Based on the expert feedback, the parameters are adjusted for the next phase, narrowing the search space.

Steps 1 and 2 are repeated until satisfactory segmentation results are achieved. In our work, we use 4 phases, i.e. 3 iterations after the first phase.
This HITL approach is designed to improve the accuracy and robustness of GFP stack segmentation by leveraging knowledge from human experts.

## Usage

1.  **Prepare Image Stacks:** Place your GFP image stacks in the `../Results/../Stacks_hotl/` directory. Organize the stacks into subdirectories corresponding to the channel names (e.g., `../Results/../Stacks_hotl/ST6GAL1/`).

2.  **Configure Parameters:** Modify the parameter dictionaries (`meanpixels`, `cutoffs`, `paramtypes`, etc.) within the script to reflect your specific image data and segmentation requirements.

3.  **Set the Phase:** Adjust the `phase` variable to indicate the current iteration of the HITL process. The fractional deviation for parameter generation is automatically calculated based on the phase.

4.  **Run the Script:** Navigate to the hitl folder and run the following script.

    ```bash
    <path-to-python>/python hitl_generation.py
    ```

5.  **Review Segmentation Results:** The segmented images will be saved in the `../Results/../<channelname>/` directory, organized by channel.

6.  **Provide Expert Feedback:** Review the segmentation results and determine if further parameter refinement is needed. Based on the review, adjust the `meanpixels` and `cutoffs` dictionaries for the next phase.

7.  **Iterate:** Repeat steps 3-6 for each phase of the HITL process.

## Configuration

The script uses the following configuration parameters, now accessible via command-line arguments:

-   **`--phase`:** Indicates the current phase of the HITL process. Start with 1.
-   **`--input_dir`:** Specifies the directory containing the original image stacks.
-   **`--output_dir`:** Specifies the directory where segmented images are saved. 
-   **`--do_channels`:** A list of channels to process.

In addition to the command-line arguments, the following parameters are still configured within the script:

-   **`meanpixels`:** A dictionary mapping channel names to mean pixel dimensions. These dimensions are used to generate scale parameters.
-   **`cutoffs`:** A dictionary mapping channel names to lists of cutoff values.
-   **`channels`:** A dictionary mapping directory names to channel names.
-   **`usefunction`:** A dictionary mapping directory names to segmentation functions.
-   **`paramtypes`:** A dictionary mapping channel names to relevant parameter types (e.g., `f2params`, `s2params`).
-   **`stgal_topothin_pvals`:** A list of parameter values for the topological thinning in ST6GAL1 channel.
