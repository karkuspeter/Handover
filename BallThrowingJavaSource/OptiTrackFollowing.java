package application;

import java.io.IOException;
import java.net.DatagramPacket;
import java.net.DatagramSocket;
import java.net.InetAddress;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import static com.kuka.roboticsAPI.motionModel.BasicMotions.ptp;

import com.kuka.common.ThreadUtil;
import com.kuka.roboticsAPI.applicationModel.RoboticsAPIApplication;
import com.kuka.roboticsAPI.controllerModel.Controller;
import com.kuka.roboticsAPI.deviceModel.LBR;
import com.kuka.roboticsAPI.geometricModel.CartDOF;
import com.kuka.roboticsAPI.geometricModel.Frame;
import com.kuka.roboticsAPI.geometricModel.ObjectFrame;
import com.kuka.roboticsAPI.geometricModel.PhysicalObject;
import com.kuka.roboticsAPI.geometricModel.Tool;
import com.kuka.roboticsAPI.geometricModel.Workpiece;
import com.kuka.roboticsAPI.geometricModel.World;
import com.kuka.roboticsAPI.motionModel.ISmartServoRuntime;
import com.kuka.roboticsAPI.motionModel.LIN;
import com.kuka.roboticsAPI.motionModel.ServoMotion;
import com.kuka.roboticsAPI.motionModel.SmartServo;
import com.kuka.roboticsAPI.motionModel.controlModeModel.CartesianImpedanceControlMode;
import com.kuka.roboticsAPI.sensorModel.ForceSensorData;
import com.kuka.roboticsAPI.userInterface.ServoMotionUtilities;
import com.kuka.generated.ioAccess.*;


public class OptiTrackFollowing extends RoboticsAPIApplication {
	private Controller kuka_Sunrise_Cabinet_1;
	private LBR lbr_iiwa_14_R820_1;
	private Tool ToolGripper;
	
	private Workpiece Object_05;

	private CartesianImpedanceControlMode mode = null ; 
	private CartesianImpedanceControlMode modeHard = null ; 
	private CartesianImpedanceControlMode modeHandShake = null ; 

	private GripperIOGroup FingerTip;
	private ATIIOGroup ForceTorqueSensor ;
	
    private UdpClient client = null;

	
	Boolean _IsRunning = true;
	double[] OptiTrackData;
	int count = 0;
	
	double robotBaseX = 0.14;
	double robotBaseY = 0.2;
	double robotBaseZ = -0.02;
			
	

	
	int timeframe ; 
	double RB1_x = 0.0 ;
	double RB1_y = 0.0 ;
	double RB1_z = 0.0 ;
	double RB1_qx = 0.0 ;
	double RB1_qy = 0.0 ;
	double RB1_qz = 0.0 ;
	double RB1_qw = 0.0 ;


	double RB2_x = 0.0 ;
	double RB2_y = 0.0 ;
	double RB2_z = 0.0 ;
	double RB2_qx = 0.0 ;
	double RB2_qy = 0.0 ;
	double RB2_qz = 0.0 ;
	double RB2_qw = 0.0 ;


	
	public enum RoboState{
		Initial,
		Approach,
		Catching,
		Delivery,
		Drop,
		ExceptionApproach,
		Final;
	}
	
	
  RoboState State = RoboState.Initial;
	
	byte[] receiveDataS = new byte[1000];

	DatagramPacket PacketsS = new DatagramPacket(receiveDataS, receiveDataS.length);

	byte[] receiveData = new byte[1024];

	DatagramPacket Packets = new DatagramPacket(receiveData, receiveData.length);
	
	//smartServo
	//private LBR _theLbr;
    //private PhysicalObject _toolAttachedToLBR;
     
    
    @Override
    public void initialize()
    {
        // Locate the "first" Lightweight Robot in the system
    	lbr_iiwa_14_R820_1 = ServoMotionUtilities.locateLBR(getContext());
        // FIXME: Set proper Weights or use the plugin feature
        
    	kuka_Sunrise_Cabinet_1 = (Controller) getContext().getControllers().toArray()[0];
//    	lbr_iiwa_14_R820_1 = (LBR) kuka_Sunrise_Cabinet_1.getDevices().toArray()[0];
        
		FingerTip= new GripperIOGroup(kuka_Sunrise_Cabinet_1);
		ForceTorqueSensor = new ATIIOGroup(kuka_Sunrise_Cabinet_1);

		
//		Gripper=getApplicationData().createFromTemplate("Gripper");
		ToolGripper=getApplicationData().createFromTemplate("ToolGripper");		
			
		mode = new CartesianImpedanceControlMode() ; 
		mode.parametrize(CartDOF.A).setStiffness(100) ;  	// rotation about Z
		mode.parametrize(CartDOF.B).setStiffness(100);		// rotation about Y
		mode.parametrize(CartDOF.C).setStiffness(100);		// rotation about X
		
		//mode.parametrize(CartDOF.TRANSL).setStiffness(4000.0) ;
		mode.parametrize(CartDOF.X).setStiffness(1000) ;
		mode.parametrize(CartDOF.Y).setStiffness(1000) ;
		mode.parametrize(CartDOF.Z).setStiffness(1000) ;
		mode.parametrize(CartDOF.ALL).setDamping(0.7) ;
		
		
		modeHandShake = new CartesianImpedanceControlMode() ; 
		modeHandShake.parametrize(CartDOF.A).setStiffness(50) ;  	// rotation about Z
		modeHandShake.parametrize(CartDOF.B).setStiffness(50);		// rotation about Y
		modeHandShake.parametrize(CartDOF.C).setStiffness(10);		// rotation about X
		
		//mode.parametrize(CartDOF.TRANSL).setStiffness(4000.0) ;
		modeHandShake.parametrize(CartDOF.X).setStiffness(100) ;
		modeHandShake.parametrize(CartDOF.Y).setStiffness(100) ;
		modeHandShake.parametrize(CartDOF.Z).setStiffness(500) ;		
		modeHandShake.parametrize(CartDOF.ALL).setDamping(0.7) ;
		
	
		
		modeHard = new CartesianImpedanceControlMode() ; 
		modeHard.parametrize(CartDOF.A).setStiffness(100) ;  	// rotation about Z
		modeHard.parametrize(CartDOF.B).setStiffness(100);		// rotation about Y
		modeHard.parametrize(CartDOF.C).setStiffness(100);		// rotation about X
		
		//mode.parametrize(CartDOF.TRANSL).setStiffness(4000.0) ;
		modeHard.parametrize(CartDOF.X).setStiffness(1000) ;
		modeHard.parametrize(CartDOF.Y).setStiffness(1000) ;
		modeHard.parametrize(CartDOF.Z).setStiffness(1000) ;
		
		modeHard.parametrize(CartDOF.ALL).setDamping(0.7) ;
		
    }
	
    
	public byte[] stringToBytes(String ss)
    {
    	ByteBuffer buff = ByteBuffer.allocate(ss.length()*2);
        buff.order(ByteOrder.LITTLE_ENDIAN);
        buff.put(ss.getBytes());
        return buff.array();
    }
	
		
	public int ReceiveData()
	{
		
		String recvStr = client.getString();
		String[] optiTrackStr = recvStr.split(" ");
		OptiTrackData=new double[optiTrackStr.length];
		for(int j=0;j<optiTrackStr.length;j++)
			{
			OptiTrackData[j]=Double.valueOf(optiTrackStr[j]);

			}
    	return 1;
	}
	
	
	
	@Override 
	public void dispose(){
		
		_IsRunning = false ; 
		System.out.println(" closing sockets in Dispose Block"); 
		
		if(client != null)
		{
			client.kill();
		}
		
		super.dispose();
	}
	
	
	/**
     * Move to an initial Position WARNING: MAKE SHURE, THAT the pose is
     * collision free.
     */
    private void moveToInitialPosition()
    {
    	lbr_iiwa_14_R820_1.move(ptp( Math.toRadians(0.0), Math.toRadians(20.0), Math.toRadians(5.0), Math.toRadians(-96.0) , Math.toRadians(2.0), Math.toRadians(-13.0), Math.toRadians(-35.0)).setJointVelocityRel(0.1));
        
    	if (FingerTip.getActReq() != 9)
    		{
    		FingerTip.setActReq(9);
    		}
    	
    	FingerTip.setSpeed(255);
		FingerTip.setForce(5);
		
        FingerTip.setPosReq(255); // close the gripper
    	
        ThreadUtil.milliSleep(2000); // wait for 2 seconds
        
    }
	
        
    private double[] getForce(double[] InitForces) {

		ForceSensorData torData = lbr_iiwa_14_R820_1.getExternalForceTorque(ToolGripper.getFrame("/BasePad"));

		double ForceX = torData.getForce().getX() - InitForces[0];
		double ForceY = torData.getForce().getY() - InitForces[1];
		double ForceZ = torData.getForce().getZ() - InitForces[2];
		
		
		double ForceMagnitude = Math.sqrt(	Math.pow(ForceX, 2) + 
											Math.pow(ForceY, 2) + 
											Math.pow(ForceZ, 2)		);

		
		double[] Forces = new double[3];
		Forces[0] = ForceX;
		Forces[1] = ForceY;
		Forces[2] = ForceZ;

		return Forces;

	}

    

    

public void run() {	
	//  RoboState State = RoboState.Initial;

		try{

		  	try {
		  		client = new UdpClient();
				client.start();

		  		ReceiveData();

		  		if(OptiTrackData[0]==-1)
		      	{
		      		System.err.println("No optitrack data");
		      		return;
		      	}
		      	else
		      	{
		      		String S_reply = "Connected";
		      		System.out.println(S_reply);
		      	}
		  	    
		      	} 
		  	
		  		catch (Exception ex) {
		        System.err.println("Can't accept client connection. ");
		  		}	
		  	
		  	
		  	ToolGripper.attachTo(lbr_iiwa_14_R820_1.getFlange());

		  	moveToInitialPosition();
		  	
	        SmartServo aSmartServoMotion = new SmartServo(lbr_iiwa_14_R820_1.getCurrentJointPosition());

	        aSmartServoMotion.useTrace(true);

	        // for Automatic mode 0.25, for T1 mode 1
	        aSmartServoMotion.setJointAccelerationRel(0.2);
	        aSmartServoMotion.setJointVelocityRel(0.2);

	        aSmartServoMotion.setMinimumTrajectoryExecutionTime(5e-3);

	        System.out.println("Starting RealtimeMotion in Position Mode");
	        
  			if (aSmartServoMotion.validateForImpedanceMode(ToolGripper) != true)
  			{
  			    getLogger().info("Validation for SmartServo Compliant control failed");
  			}
	        
	        
	        //_toolAttachedToLBR.getDefaultMotionFrame().moveAsync(aSmartServoMotion.setMode(controlMode));
	        
	        ToolGripper.getFrame("/BasePad").moveAsync(aSmartServoMotion.setMode(modeHandShake));

	        // Fetch the Runtime of the Motion part
	        ISmartServoRuntime theServoRuntime = aSmartServoMotion.getRuntime();

	        //Frame aFrame = theServoRuntime.getCurrentCartesianDestination(Gripper.getDefaultMotionFrame());

		      //. Receive data and move the robot
		  	while(_IsRunning)
		  	{
		  		count = count + 1;
		  		
		  		ReceiveData();
		  		
		  		
		  		if(OptiTrackData[0]==-1)
			  	{
		  			System.out.println(" stop motion ! at " + theServoRuntime.toString());
		  	        // Stop the motion
		  	        theServoRuntime.stopMotion();
		  	        
			  		break;
			  	}
		  		
		  		
		  		if(OptiTrackData[0]!=-1)
		  		{
		  		
			  			//smartServo
				  	    // Synchronize with the realtime system
				  		theServoRuntime.updateWithRealtimeSystem();
		  			
			  			
			  			// OptiTrackData : 1st element is time frame stamp
				  		// followed by rigidBody data : for each rigidBody: RB_index , x , y , z , qx , qy , qz , qw 
			  
				  		timeframe = (int)OptiTrackData[0];
				  		
				  		//RB1_index = (int) OptiTrackData[1];
			  			RB1_x = (OptiTrackData[2] - robotBaseX)*1000;
			  			RB1_y = (OptiTrackData[3] - robotBaseY)*1000;
			  			RB1_z = (OptiTrackData[4] - robotBaseZ)*1000;
			  			RB1_qx = OptiTrackData[5];
			  			RB1_qy = OptiTrackData[6];
			  			RB1_qz = OptiTrackData[7];
			  			RB1_qw = OptiTrackData[8];
			  			
				  		//RB2_index = (int) OptiTrackData[9];
			  			RB2_x = (OptiTrackData[10] - robotBaseX)*1000;
			  			RB2_y = (OptiTrackData[11] - robotBaseY)*1000;
			  			RB2_z = (OptiTrackData[12] - robotBaseZ)*1000;
			  			RB2_qx = OptiTrackData[13];
			  			RB2_qy = OptiTrackData[14];
			  			RB2_qz = OptiTrackData[15];
			  			RB2_qw = OptiTrackData[16];
			  			
			  			
			  			//System.out.println("RB1 XYZ = " + String.valueOf( RB1_x ) + ", "+ String.valueOf( RB1_y ) + ", "+ String.valueOf( RB1_z ));
		            			            	
		            	Frame msrPoseUpdate = new Frame();
		            	msrPoseUpdate = theServoRuntime.getCurrentCartesianPosition(ToolGripper.getFrame("/BasePad/FingerTipStraight"));
		            	
		            	double posXmsr = msrPoseUpdate.getX();
		            	double posYmsr = msrPoseUpdate.getY();
		            	double posZmsr = msrPoseUpdate.getZ();
		            	
		           		double oriAlpha = msrPoseUpdate.getAlphaRad();   	// rotation about Z
        		  		double oriBeta  = msrPoseUpdate.getBetaRad();		// rotation about Y
        		  		double oriGamma = msrPoseUpdate.getGammaRad();		// rotation about X
		             	
        		  		// read the force applied on the gripper
        		  		double[] InitForces = new double[3];
        				InitForces[0]= 0;
        				InitForces[1]= 0;
        				InitForces[2]= 0;
        				
        				double[] ForceOnGripper = getForce(InitForces);
		            	double Fx = ForceOnGripper[0];
		            	double Fy = ForceOnGripper[1];
		            	double Fz = ForceOnGripper[2];

		            	//System.out.println("Fx= " + String.valueOf(Fx) + ", Fy= " + String.valueOf(Fy) + ", Fz= " + String.valueOf(Fz)  );	
		            	
			  			
		            	switch(State){
	                	
			  			case Initial:
	                		
			  				Frame InitFrame =  lbr_iiwa_14_R820_1.getCurrentCartesianPosition(ToolGripper.getFrame("/BasePad/FingerTipStraight"));
	                		
			  	          
	                		InitFrame.setX(350);
	                    	InitFrame.setY(200);
	                    	InitFrame.setZ(400);
	                    	
//	                    	InitFrame.setAlphaRad( Math.toRadians(-90) 	);
//	                    	InitFrame.setBetaRad(  Math.toRadians(0) 	);
//	                    	InitFrame.setGammaRad( Math.toRadians(180) 	);
	                    	
	                    	
	                    	FingerTip.setPosReq(200);

//	                    	theServoRuntime.setDestination(InitFrame,World.Current.getRootFrame());
	                    	ThreadUtil.milliSleep(1000);
	                    	// go to Drop position

			  				State = RoboState.Approach;
	                		System.out.println("Init state done => Approach state");
	                		System.out.println("goal: x = " + String.valueOf(RB1_x) + " , y = " + String.valueOf(RB1_y) + " , z = " + String.valueOf(RB1_z) );

	                		break;
	                	
	                	
			  			
			  			case Approach:
                		
			  				try {
                    	
        		  
	        		  		Frame goalFrame = new Frame();  
	        		  		
	        				goalFrame = lbr_iiwa_14_R820_1.getCurrentCartesianPosition(ToolGripper.getFrame("/BasePad/FingerTipStraight"));
                		
        				

        					double Xlimit = 1500.0 ; 
        					double Ylimit = 1500.0 ; 
        					double ZlimitLow = -100.0 ; 
        					double ZlimitHigh = 800 ; 

        				
                            // boxing bound, saturate the set points (RB) to the limits  
                            if ( Math.abs(RB1_x) > Xlimit) 
                            {
                            	RB1_x = Xlimit * Math.signum(RB1_x);
                            }
                            
                            if ( Math.abs(RB1_y) > Ylimit)
                            {
                            	RB1_y = Ylimit * Math.signum(RB1_y);
                            }
                            
                            // avoid hitting the table
                            if ( Math.abs(RB1_z) <= ZlimitLow) 
                            {
                            	RB1_z = ZlimitLow  ; 
                            }
                            
                            // 
                            if ( Math.abs(RB1_z) >= ZlimitHigh) 
                            {
                            	RB1_z = ZlimitHigh  ; 
                            }
                                                        
                            // compute the absolute motion w.r.t. World 
            		  		goalFrame.setX(RB1_x);
                            goalFrame.setY(RB1_y);
                            goalFrame.setZ(RB1_z);
//                            goalFrame.setAlphaRad(oriAlpha);
//                            goalFrame.setBetaRad(oriBeta);
//                            goalFrame.setGammaRad(oriGamma);
                            
//                            goalFrame.setAlphaRad(-90 * Math.PI / 180);
//                            goalFrame.setBetaRad (  0 * Math.PI / 180);
//                            goalFrame.setGammaRad(180 * Math.PI / 180);
                            theServoRuntime.setDestination(goalFrame, World.Current.getRootFrame());

                            double distance = Math.sqrt(  		Math.pow( (posXmsr - RB1_x), 2) 
						    								+ 	Math.pow( (posYmsr - RB1_y), 2) 
															+ 	Math.pow( (posZmsr - RB1_z), 2) 	) ; 
                            
                            
//                            double distance_ball = Math.sqrt( Math.pow( (posXmsr - OptiTrackData[11]), 2) 
//    								+ Math.pow( (posYmsr - OptiTrackData[12]), 2) 
//									+ Math.pow( (posZmsr - OptiTrackData[13]), 2) 	) ; 
                            
                            
                            //System.out.println("Dist=" + String.valueOf(distance));

			  				}
			  			 	
	            		  	catch (Exception ex2) 
	            		          {
	            		    	State = RoboState.ExceptionApproach;
	            		     
	            		    	System.out.println("ExceptionApproach ");
	            		  		
	            		          }
			  				
			  				break;
	
                	
                	case Final:	
                		
                
                		
                		
                		ThreadUtil.milliSleep(2500);
                		State = RoboState.Initial;
    
                		System.out.println("Initial state");
                		
                		break;
                		
                		
                		
                        case ExceptionApproach:	
                		
                        ThreadUtil.milliSleep(100);
                		
                		State = RoboState.Approach;
    
//                		System.out.println("Approach state");
                		
                		break;
                	}     //switch
		  			
		    
                    } 
                    
                    
            }
		  		
		  		
		  	
		  	

		}
		catch(Exception e){
		      e.printStackTrace();
		  }

	
	}
	

	
	/**
	 * Auto-generated method stub. Do not modify the contents of this method.
	 */
	public static void main(String[] args) {
		OptiTrackFollowing app = new OptiTrackFollowing();
		app.runApplication();
	}
}
