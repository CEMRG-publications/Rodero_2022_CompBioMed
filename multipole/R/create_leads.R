source('./preprocessing.R')

for(heart_num in c(1:24)){
	    if(heart_num <= 20){
		            Create_n_electrodes(heart_num,"h","AN")
        Create_n_electrodes(heart_num,"h","AL")
	        Create_n_electrodes(heart_num,"h","LA")
	        Create_n_electrodes(heart_num,"h","IL")
		        Create_n_electrodes(heart_num,"h","IN")
		    }

    Create_n_electrodes(heart_num,"HF","AN")
        Create_n_electrodes(heart_num,"HF","AL")
        Create_n_electrodes(heart_num,"HF","LA")
	    Create_n_electrodes(heart_num,"HF","IL")
	    Create_n_electrodes(heart_num,"HF","IN")
}

system("/home/crg17/Desktop/scripts/4chmodel/sh/sendmail.sh \"All leads created\"")
