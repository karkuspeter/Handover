package application;




import java.io.IOException;
import java.net.DatagramPacket;
import java.net.DatagramSocket;
import java.net.InetAddress;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.concurrent.TimeUnit;

import com.kuka.common.ThreadUtil;
import com.kuka.roboticsAPI.applicationModel.RoboticsAPIApplication;
import com.kuka.roboticsAPI.controllerModel.Controller;
import com.kuka.roboticsAPI.deviceModel.JointPosition;
import com.kuka.roboticsAPI.deviceModel.LBR;
import com.kuka.roboticsAPI.geometricModel.Frame;
import com.kuka.roboticsAPI.geometricModel.ObjectFrame;
import com.kuka.roboticsAPI.geometricModel.PhysicalObject;
import com.kuka.roboticsAPI.geometricModel.Tool;
import com.kuka.roboticsAPI.geometricModel.Workpiece;

import com.kuka.roboticsAPI.motionModel.ISmartServoRuntime;
import com.kuka.roboticsAPI.motionModel.LIN;
import com.kuka.roboticsAPI.motionModel.SmartServo;
import com.kuka.roboticsAPI.motionModel.SplineMotionCP;
import com.kuka.roboticsAPI.motionModel.controlModeModel.CartesianImpedanceControlMode;
import com.kuka.roboticsAPI.userInterface.ServoMotionUtilities;
import com.kuka.generated.ioAccess.*;


import static com.kuka.roboticsAPI.motionModel.BasicMotions.*;


public class Masood_udp extends RoboticsAPIApplication {
	private Controller kuka_Sunrise_Cabinet_1;
	private Tool ToolGripper;
	private Workpiece Pen;
	private CartesianImpedanceControlMode mode = null ; 
	private GripperIOGroup FingerTip;
	private LBR lbr_iiwa_14_R820_1;
	
	Boolean _IsRunning = true;
	byte[] value_frame;
	
	double posX_ini = 0;
	double posY_ini = 0;
	double posZ_ini = 0;
	char gripperCmd = 'c';//Command for gripper
	char gripperSts = 'o';//Gripper status
	byte[] receiveDataS = new byte[5];
	double[] gripper_hand = new double[2]; // Distance from gripper to left & right hands
   
    private JointPosition MoveToPosition;
    private UdpClient client = null;
    
    public void initialize()
    {
        // Locate the "first" Lightweight Robot in the system
//        _theLbr = ServoMotionUtilities.locateLBR(getContext());
        // FIXME: Set proper Weights or use the plugin feature
//        lbr_iiwa_14_R820_1 = ServoMotionUtilities.locateLBR(getContext());
    	kuka_Sunrise_Cabinet_1 = (Controller) getContext().getControllers().toArray()[0];
    	lbr_iiwa_14_R820_1 = (LBR) kuka_Sunrise_Cabinet_1.getDevices().toArray()[0];
		FingerTip= new GripperIOGroup(kuka_Sunrise_Cabinet_1);
		ToolGripper = getApplicationData().createFromTemplate("ToolGripper");
//		Pen=getApplicationData().createFromTemplate("Pen");
		System.out.println("initialize end!" );
//		lbr_iiwa_14_R820_1.move(lin(640.0,325.0,290.0));

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
	

	
	public void run() {	
		

		try{
	
			client = new UdpClient();
			client.start();
			
		  	FingerTip.setActReq(9);
			FingerTip.setSpeed(120);
			FingerTip.setForce(10);
			FingerTip.setPosReq(0) ; //Open
			ToolGripper.attachTo(lbr_iiwa_14_R820_1.getFlange());

			System.out.println("Usage: close hand when object 1 pos x < 2.0m ");
			System.out.println("open hand when object 1 pos x >= 2.0m ");
			System.out.println("exit when object 1 pos x >= 3.0m ");

//        	theServoRuntime.setDestination(destFrame, World.Current.getRootFrame());

			while(_IsRunning)
			{

				Frame destination = lbr_iiwa_14_R820_1.getCurrentCartesianPosition(ToolGripper.getFrame("/BasePad"));

				String recvStr = client.getString();

				//System.out.printf("l: %d\n" , recvPacket.getLength() );
				//System.out.println("receive data:" + recvStr);

				if (recvStr.length() > 0){
					String[] optiTrackStr = recvStr.split(" ");
					double[] optiTrackDouble=new double[optiTrackStr.length];
					for(int j=0;j<optiTrackStr.length;j++)
					{
						optiTrackDouble[j]=Double.valueOf(optiTrackStr[j]);
						//System.out.printf("float output:%f\n" , optiTrackDouble[j]);
					}
					if (optiTrackDouble[2] > 3.0){ //pos x for object 1
						break;
					}
					if (optiTrackDouble[2] < 2.0) {
						FingerTip.setPosReq(255) ; //close
					} else
					{
						FingerTip.setPosReq(0) ; //open
					}
				}


				//ThreadUtil.milliSleep(500);

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
		Masood_udp app = new Masood_udp();
		app.runApplication();
	}
}