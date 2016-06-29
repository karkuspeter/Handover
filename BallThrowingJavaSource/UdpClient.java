package application;

import java.net.*;
import java.util.*;
import java.io.*;

public class UdpClient extends Thread {
	protected DatagramSocket socket = null;
	protected String recvStr = new String();
	protected Boolean end = false;
	protected long LastAccess = System.currentTimeMillis();
	protected List<String> buffer;
	protected Boolean isRecording = false;
	
		
	public String getString(){
		String s;
		synchronized(this){
			s = recvStr;
			LastAccess = System.currentTimeMillis();
		}
		return s;
	}
	
	public UdpClient() throws IOException{
		this("UdpClient");
	}
	
	public UdpClient(String name) throws IOException{
		super(name);
	  	
	}
	public void kill()
	{
		end = true;
	}

	public void dispose() {
		if(socket != null)
		{	
			if(!(socket.isClosed()))
			{
				socket.close();
			}
		}

		System.out.println("Disposing Client");
	}

	public void StartRecording(){
		synchronized(this){
			buffer.clear();
			isRecording = true;
		}
		
	}
	
	public List<String> GetRecording(){
		synchronized(this){
			isRecording = false;
		}
		return buffer;
	}
	
	public void run() {
		   byte[] recvBuf = new byte[2000];
		   DatagramPacket recvPacket;
		  	try {
		  		socket = new DatagramSocket(12358);
		    } catch (IOException ex) {
		          System.err.println("Can't setup server on this port number. ");
		          System.err.println(ex);
		    }		
		  	System.out.printf("getLocalPort:%d\n" ,socket.getLocalPort() );
		  	System.out.printf("getLocalAddress:%s\n" ,socket.getLocalAddress() );
		   
		   while( !end ){
			   recvPacket = new DatagramPacket(recvBuf , recvBuf.length);
				try{
					socket.receive(recvPacket);
				}
				catch(Exception e){
					System.out.println("socket close error. ");
					e.printStackTrace();
					break;
				}

				synchronized(this){
					recvStr = new String(recvPacket.getData() , 0 , recvPacket.getLength());
					if (isRecording) {
						buffer.add(recvStr);
					}
					if (System.currentTimeMillis()-LastAccess > 10000){
						break;
					}
				}
		   }
		   
			if(socket != null)
			{	
				if(!(socket.isClosed()))
				{
					socket.close();
				}
			}		   
		   System.out.println("Exiting Client");

	}
}