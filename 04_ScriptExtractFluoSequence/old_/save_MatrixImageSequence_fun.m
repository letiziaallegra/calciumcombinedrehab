function ImageSequence = save_MatrixImageSequence_fun(Name,Date,MatrixImageSequence,CentralFrame,ROISelected,BeforeCentralFrame,AfterCentralFrame,downsamplingfactor)


ImageSequence.Name                = Name;
ImageSequence.Date                = Date;
ImageSequence.MatrixImageSequence = MatrixImageSequence;
ImageSequence.CentralFrameInfo    = CentralFrame;
ImageSequence.ROISelected         = ROISelected;
ImageSequence.CentFrame           = BeforeCentralFrame+1;
ImageSequence.BeforeCentralFrame  = BeforeCentralFrame;
ImageSequence.AfterCentralFrame   = AfterCentralFrame;
ImageSequence.downsamplingfactor  = downsamplingfactor;


Filename = ['ImageSequence_',Name,'_',Date,'_CenFrame_',num2str(CentralFrame(1)),'_',num2str(CentralFrame(2)),'_ROI_',num2str(ROISelected)];
save(Filename,'ImageSequence','-v7.3');
                         
end