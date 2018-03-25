function ImageSequence = save_MatrixImageSequence_fun_Rest(Name,Date,MatrixImageSequence,CentralFrame,ROISelected,BeforeCentralFrame,AfterCentralFrame,downsamplingfactor)


ImageSequence.Name                = Name;
ImageSequence.Date                = Date;
ImageSequence.MatrixImageSequence = MatrixImageSequence;
ImageSequence.CentralFrameInfo    = CentralFrame;
ImageSequence.ROISelected         = ROISelected;
ImageSequence.CentFrame           = CentralFrame;
ImageSequence.BeforeCentralFrame  = BeforeCentralFrame;
ImageSequence.AfterCentralFrame   = AfterCentralFrame;
ImageSequence.downsamplingfactor  = downsamplingfactor;


Filename = ['ImageSequence_',Name,'_',Date,'_Rest','_ROI_',num2str(ROISelected)];
save(Filename,'ImageSequence','-v7.3');
                         
end