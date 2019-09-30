# EM-IQM
Implementation of "Image Quality Assessment for DIBR Synthesized Views using Elastic Metric"


# Running:
em_score = EM_IQM( im_ori,im_syth , num_fp)

- input
im_ori : the reference/original image
im_syn : the synthesized/distorted image 
num_fp : number of feature pointes, in the paper it is set as num_fp=500

-output 
em_score : the predicted quality score

# Citation
@inproceedings{ling2017image,
  title={Image quality assessment for dibr synthesized views using elastic metric},
  author={Ling, Suiyi and Le Callet, Patrick},
  booktitle={Proceedings of the 25th ACM international conference on Multimedia},
  pages={1157--1163},
  year={2017},
  organization={ACM}
}
