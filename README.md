
# Our_p3p
Codes (encrypted using MATLAB pcode()) and demos for our p3p method, the codes will be released after being acceptted.

# Thirdparty Methods (Folder "P3P_tool" )
The thirdparty p3p methods are from the p3p-toolbox provided by Gaku Nakano. Due to the license, we do not plan to provide the toolbox here.  You can download it from [link](https://github.com/g9nkn/p3p_problem) and then decompress it into the "P3P_tool" folder. After that, you will be able to run the following demos:
 - demo_test_noise_robustness
 - demo_test_numerical_stability_rt
 - demo_test_runtime

To run the "demo_test_ETH3D_VGG", you also need to add the following codes after Lines 209 (i.e., after `Qs = reshape(Qs,4,4);`) in the "p3p_gao" function of the toolbox to avoid NAN-error:

    if all(all(isnan(Qs)))
        R = nan(3);  t = nan(3,1);
        return
    end


# Demo Results (Folder "demo_results")
## Numerical Stability
![ns-R](https://github.com/Johnnyzyzy/our_p3p/blob/main/demo_results/ns-R.png)
![ns-t](https://github.com/Johnnyzyzy/our_p3p/blob/main/demo_results/ns-t.png)

## Noise Robustness
![noise-robustness-R](https://github.com/Johnnyzyzy/our_p3p/blob/main/demo_results/nr-R.png)
![noise-robustness-t](https://github.com/Johnnyzyzy/our_p3p/blob/main/demo_results/nr-t.png)

## Runtime
**Note that the runtime varies with different hard- or soft-wares or systems.**

Screen recording video: ./demo_results/rumtime.mp4

Ours experiment environments:
 - PC: Xiaomi Mi Notebook Pro 15 OLED Enhanced Edition
 - CPU: 11th Gen Intel(R) Core(TM) i7-11390H @ 3.40GHz   3.42 GHz
 - RAM: 16 GB
 - SYSTEM: WIN10 X64
 - GPU: NVIDIA®  GeForce®  MX450

![runtime](https://github.com/Johnnyzyzy/our_p3p/blob/main/demo_results/runtime.png)
![runtime-gif1](https://github.com/Johnnyzyzy/our_p3p/blob/main/demo_results/runtime1.gif)
![runtime-gif2](https://github.com/Johnnyzyzy/our_p3p/blob/main/demo_results/runtime2.gif)
