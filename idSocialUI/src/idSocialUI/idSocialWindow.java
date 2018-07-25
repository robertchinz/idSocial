package idSocialUI;

import java.awt.EventQueue;

// import idSocialOptionsEditor.RowEditorModel;
import javax.swing.*;
import java.io.*;
import java.awt.Window.Type;

import javax.swing.table.DefaultTableModel;

import java.awt.Font;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

import net.miginfocom.swing.MigLayout;

import javax.swing.GroupLayout.Alignment;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.event.TreeSelectionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import javax.swing.tree.TreePath;
import javax.swing.filechooser.*;

import java.util.List;
import java.util.ArrayList;


public class idSocialWindow extends JFrame {

	// private JFrame frmIdsocial;
	private JTable table;
	private JTable table_1;
	private IdSocialOptionTable table_trinfo;
	private IdSocialOptionTable table_troptions;
	

	private JButton btnNewButton;

	private JTree tree = new JTree();
	private JTree trajTree = new JTree();
	private Object[] TrajectoryList; // List converted from file tree containing
										// trajectory locations
	
	private File workingDirectory = new File(System.getProperty("user.dir"));

	
	public IdSocialOptionTable getTrajectoryOptions() 
	{
		return table_troptions;
	}
	
	
	public JButton getPreprocessButton()
	{
		return  btnNewButton;
	}

	public JTree getTrajectoryTree()
	{
		return  trajTree;
	}

	/**
	 * Launch the application.
	 */
	public static void main(String[] args) {

		

		try {

			// String lookAndFeelString =
			// "javax.swing.plaf.nimbus.NimbusLookAndFeel";
			// UIManager.setLookAndFeel(lookAndFeelString);
			UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
		}

		catch (ClassNotFoundException | InstantiationException | IllegalAccessException
				| UnsupportedLookAndFeelException e) {
			System.err.println("Couldn't find class for specified look and feel");
			System.err.println("Did you include the L&F library in the class path?");
			System.err.println("Using the default look and feel.");
		}

		EventQueue.invokeLater(new Runnable() {
			public void run() {

				try {
					idSocialWindow window = new idSocialWindow();
					window.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	/**
	 * Create the application.
	 */
	public idSocialWindow() {
		initialize();
	}

	/**
	 * Initialize the contents of the frame.
	 */
	private void initialize() {

		// frmIdsocial = new JFrame();
		setFont(new Font("Dialog", Font.PLAIN, 14));
		setTitle("idSocial");
		setBounds(100, 100, 968, 569);
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);

		JMenu mnFile = new JMenu("File");
		menuBar.add(mnFile);

		JMenu mnView = new JMenu("View");
		menuBar.add(mnView);
		getContentPane().setLayout(
				new MigLayout("", "[150px,grow][290px,grow][112px]", "[337px,grow][323px,grow]"));

		// File tree for groups, days, trials .....
		JScrollPane scrollPane = new JScrollPane();
		getContentPane().add(scrollPane, "cell 0 0 1 1,grow");

		DefaultMutableTreeNode rootNode = new DefaultMutableTreeNode("Project");
		DefaultMutableTreeNode group1 = new DefaultMutableTreeNode("Add Group");

		rootNode.add(group1);

		trajTree.setModel(new DefaultTreeModel(rootNode));
		scrollPane.setViewportView(trajTree);

		trajTree.setRootVisible(false);
		trajTree.setEditable(true);

		trajTree.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent me) {
				trajTreeClick(me);
			}
		});

		// Info table
//		JScrollPane scrollPane_trinfo = new JScrollPane();
//
//		getContentPane().add(scrollPane_trinfo, "cell 1 0,grow");
//
//		table_trinfo = new IdSocialOptionTable();
//		scrollPane_trinfo.setViewportView(table_trinfo);
//		table_trinfo.setModel(new DefaultTableModel(new Object[][] {}, new String[] { "Information" }));
//		table_trinfo.getColumnModel().getColumn(0).setResizable(false);
		// End Plot mode table

		// Options table

		
//
//		Object[] p = {"1","2"};
//		Object[][] parValues = new Object[][] {{1,1},{true,true}};
//	    table_troptions = new IdSocialOptionTable(parValues,p);
	    table_troptions = new IdSocialOptionTable();

//	    table_troptions.setTableData(new String[] {"test", "test2","test3","test4","test5"}, new Object[] {true,1.1,new String[] {"1","2"}, new double[] {1, 2}, new double[][] {new double[]{1,2},new double[]{3,4}}});
	    table_troptions.setTableData(new String[] {"test", "test2","test3","test4"}, new Object[] {true,1.1,new String[] {"1","2"}, new double[] {1, 2}});

//	    table_troptions.getTableData() ;
//	    table_troptions = new IdSocialOptionTable(); 

//	    
	    JScrollPane scrollPane_troptions = new JScrollPane(table_troptions);

	    getContentPane().add(scrollPane_troptions, "cell 1 0 2 1,grow");

		// End Plot mode table

		btnNewButton = new JButton("Preprocess");
		btnNewButton.setFont(btnNewButton.getFont().deriveFont(btnNewButton.getFont().getSize() + 5f));
		btnNewButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
//				System.out.println("Preprocess");
//				getTrajectoryList();
				System.out.println(table_troptions.getTableData());
//				getTrajectoryList();
			}
		});
		
		
		getContentPane().add(btnNewButton, "cell 3 0,growx,aligny top");

		JScrollPane scrollPane_2 = new JScrollPane();
		DefaultMutableTreeNode methods = new DefaultMutableTreeNode("Methods");

		tree.setModel(new DefaultTreeModel(methods));

		// JTree tree = new JTree();

		scrollPane_2.setViewportView(tree);
		tree.setFont(tree.getFont().deriveFont(tree.getFont().getSize() + 5f));
		getContentPane().add(scrollPane_2, "cell 0 1,grow");

		// Plot mode table
		JScrollPane scrollPane_1 = new JScrollPane();
		table_1 = new JTable();
		scrollPane_1.setViewportView(table_1);
		table_1.setModel(new DefaultTableModel(new Object[][] {}, new String[] { "1" }));
		getContentPane().add(scrollPane_1, "cell 1 1 2 1,grow");
		table_1.getColumnModel().getColumn(0).setResizable(false);
		// End Plot mode table
		
//		String[] p = {"1","2"};
//		Object[] parValues = {1,true};
//		setTrajectoryOptionsData(p, parValues); 

		pack();
	}

	public void setTreeData(String[] parNames) {

		DefaultMutableTreeNode methods = new DefaultMutableTreeNode("Methods");
		for (int actRow = 0; actRow <= parNames.length - 1; actRow++) {

			methods.add(new DefaultMutableTreeNode(parNames[actRow]));

		}
		tree.setModel(new DefaultTreeModel(methods));
	}
	
//	public void setTrajectoryOptionsData(String[] parNames, Object[] parValues) {
//		DefaultTableModel model = (DefaultTableModel) table_troptions.getModel();
//		
//		model.setColumnIdentifiers(new Object[] { "Parameter", "Value" });
//		table_troptions.setModel(model);
//		// TableModel model = table.gemodelodel();
//		RowEditorModel rm = table_troptions.getRowEditorModel();
//
//		for (int actRow = 0; actRow <= parNames.length - 1; actRow++) {
//			
//			if (model.getRowCount()<=actRow){
//				model.addRow(new Object[]{"", ""});
//			}
//			model.setValueAt(parNames[actRow], actRow, 0);
//
//			Object elem = parValues[actRow];
//			if (elem != null) {
//				Class actClass = elem.getClass();
//				if (elem instanceof String[]) {
//					String[] combolist = (String[]) elem;
//					model.setValueAt(combolist[0], actRow, 1);
//
//					JComboBox cb = new JComboBox(combolist);
//					// cb.setSelectedItem(1);
//					((JLabel) cb.getRenderer()).setHorizontalAlignment(SwingConstants.CENTER);
//					DefaultCellEditor ed = new DefaultCellEditor(cb);
//					// inform the RowEditorMode of the situation
//					rm.addEditorForRow(actRow, ed);
//
//				} else if (elem instanceof Boolean) {
//					model.setValueAt(((Boolean) elem).booleanValue(), actRow, 1);
//					JCheckBox chkBox = new JCheckBox("", ((Boolean) elem).booleanValue());
//					chkBox.setHorizontalAlignment(JCheckBox.CENTER);
//					DefaultCellEditor ed = new DefaultCellEditor(chkBox);
//					// inform the RowEditorMode of the situation
//					rm.addEditorForRow(actRow, ed);
//
//				} else {
//					model.setValueAt(parValues[actRow], actRow, 1);
//					System.out.println(elem.getClass());
//
//				}
//				
//
//			}
//		}
//		
//	}

	public void trajTreeClick(MouseEvent me) {
		TreePath tp = trajTree.getPathForLocation(me.getX(), me.getY());
		if (tp != null) {
			// DefaultMutableTreeNode groupNew = new
			// DefaultMutableTreeNode("Group 1");
			DefaultTreeModel model = (DefaultTreeModel) trajTree.getModel();
			DefaultMutableTreeNode rootNode = (DefaultMutableTreeNode) model.getRoot();
			DefaultMutableTreeNode groupNew = new DefaultMutableTreeNode();

			// System.out.println("Selected: "+tp.getLastPathComponent());
			//
			// System.out.println("Path(1): "+tp.getPathComponent(1));

			// Adding groups
			if ((tp != null) & (tp.getPathCount() == 2) & tp.getLastPathComponent().toString().equals("Add Group")) {
				trajTree.setEditable(false);
				DefaultMutableTreeNode dayNew = new DefaultMutableTreeNode("Add Day");
				int noGroups = trajTree.getModel().getChildCount(rootNode);

				int maxGroup = 1;
				int actIndex;
				for (int actGroup = 0; actGroup < noGroups - 1; actGroup++) {
					DefaultMutableTreeNode actChild = (DefaultMutableTreeNode) rootNode.getChildAt(actGroup);
					System.out.println(actGroup);
					System.out.println(actChild);
					actIndex = buildNumber(actChild.toString());
					if (actIndex > 0) {
						maxGroup = Math.max(actIndex, maxGroup) + 1;
					} else {
						maxGroup = -1;
					}
				}

				if (maxGroup > 0) {
					groupNew = new DefaultMutableTreeNode("Group " + maxGroup);
				} else {
					groupNew = new DefaultMutableTreeNode("New Group");
				}

				groupNew.add(dayNew);
				rootNode.add(groupNew);
				model.insertNodeInto(groupNew, rootNode, noGroups - 1);
				TreePath newPath = new TreePath(groupNew.getPath());
				trajTree.expandPath(newPath);


			}
			// Adding days
			else if ((tp != null) & (tp.getPathCount() == 3) & tp.getLastPathComponent().toString().equals("Add Day")) {
				trajTree.setEditable(false);

				DefaultMutableTreeNode dayNew = new DefaultMutableTreeNode();
				DefaultMutableTreeNode groupNode = (DefaultMutableTreeNode) tp.getPathComponent(1);
				DefaultMutableTreeNode trialNew = new DefaultMutableTreeNode("Add Trial");

				int noDays = trajTree.getModel().getChildCount(groupNode);

				int maxDay = 1;
				int actIndex;

				for (int actDay = 0; actDay < noDays - 1; actDay++) {
					DefaultMutableTreeNode actChild = (DefaultMutableTreeNode) groupNode.getChildAt(actDay);
					System.out.println(actDay);
					System.out.println(actChild);
					actIndex = buildNumber(actChild.toString());
					if (actIndex > 0) {
						maxDay = Math.max(actIndex, maxDay) + 1;
					} else {
						maxDay = -1;
					}
				}
				if (maxDay > 0) {
					dayNew = new DefaultMutableTreeNode("Day " + maxDay);
				} else {
					dayNew = new DefaultMutableTreeNode("New Subgroup");
				}

				dayNew.add(trialNew);
				groupNode.add(dayNew);
				model.insertNodeInto(dayNew, groupNode, noDays - 1);
				TreePath newPath = new TreePath(dayNew.getPath());
				trajTree.expandPath(newPath);


			} else if ((tp != null) & (tp.getPathCount() == 4)
					& tp.getLastPathComponent().toString().equals("Add Trial")) {
				trajTree.setEditable(false);
				DefaultMutableTreeNode trialNew = new DefaultMutableTreeNode();
				DefaultMutableTreeNode fileNew = new DefaultMutableTreeNode();
				
				// Open file dialog
				JFileChooser fc = new JFileChooser();
				fc.setCurrentDirectory(workingDirectory);
				
				int returnVal = fc.showOpenDialog(trajTree);

				if (returnVal == JFileChooser.APPROVE_OPTION) {
					File file = fc.getSelectedFile();
					fileNew.setUserObject(file);

					workingDirectory = file.getAbsoluteFile();

				} else {

				}

				
				DefaultMutableTreeNode dayNode = (DefaultMutableTreeNode) tp.getPathComponent(2);
				int noTrials = trajTree.getModel().getChildCount(dayNode);

				int maxTrial = 1;
				int actIndex;

				for (int actTrial = 0; actTrial < noTrials - 1; actTrial++) {
					DefaultMutableTreeNode actChild = (DefaultMutableTreeNode) dayNode.getChildAt(actTrial);
					System.out.println(actTrial);
					System.out.println(actChild);
					actIndex = buildNumber(actChild.toString());
					if (actIndex > 0) {
						maxTrial = Math.max(actIndex, maxTrial) + 1;
					} else {
						maxTrial = -1;
					}
				}
				if (maxTrial > 0) {
					trialNew = new DefaultMutableTreeNode("Trial " + maxTrial);
				} else {
					trialNew = new DefaultMutableTreeNode("New Trial");
				}
				
				
				trialNew.add(fileNew);
				dayNode.add(trialNew);
				model.insertNodeInto(trialNew, dayNode, noTrials - 1);
				TreePath newPath = new TreePath(trialNew.getPath());
				trajTree.expandPath(newPath);

			} else {
				trajTree.setEditable(true);
			}

		}

	}

	private Integer buildNumber(String str) {
		StringBuilder strBuilder = new StringBuilder();
		char ch = str.charAt(str.length() - 1);

		if (Character.isDigit(ch)) {
			strBuilder.insert(0, ch);
		}

		int count = str.length() - 2;
		while (count >= 0) {
			ch = str.charAt(count);
			count = count - 1;
			if (Character.isDigit(ch)) {
				strBuilder.insert(0, ch);
			}
		}

		String digitString = strBuilder.toString();
		if (digitString.length() > 0) {
			return Integer.parseInt(digitString);
		} else {
			return -1;
		}
	}

	public String[][][][] getTrajectoryList() {
		DefaultTreeModel model = (DefaultTreeModel) trajTree.getModel();
//		List<List<String>> TrajectoryList = new List<List<String>>();
		

		DefaultMutableTreeNode rootNode = (DefaultMutableTreeNode) model.getRoot();
		int noGroups = rootNode.getChildCount();
		String[][][][] TrajectoryList = new String[noGroups-1][][][];

//		System.out.println(noGroups);
		for (int actGroup = 0; actGroup < noGroups-1; actGroup++) {
			DefaultMutableTreeNode grNode = (DefaultMutableTreeNode) rootNode.getChildAt(actGroup);
			int noDays = grNode.getChildCount();
			// System.out.println(grNode);
			TrajectoryList[actGroup] = new String[noDays-1][][];

			if (noDays > 1) {
				for (int actDay = 0; actDay < noDays-1; actDay++) {
					DefaultMutableTreeNode dayNode = (DefaultMutableTreeNode) grNode.getChildAt(actDay);
					int noTrials = dayNode.getChildCount();
					TrajectoryList[actGroup][actDay] = new String[noTrials-1][];

					// System.out.println(dayNode);

					if (noTrials>1){
					for (int actTr = 0; actTr < noTrials-1; actTr++) {
						DefaultMutableTreeNode trNode = (DefaultMutableTreeNode) dayNode.getChildAt(actTr);
						DefaultMutableTreeNode fileNode = (DefaultMutableTreeNode) trNode.getChildAt(0);

						TrajectoryList[actGroup][actDay][actTr] = new String[1];
						TreePath path = new TreePath(trNode.getPath());
						String pathString = path.toString();
//						TrajectoryList[actGroup][actDay][actTr][1] = pathString;
						TrajectoryList[actGroup][actDay][actTr][0] = fileNode.toString();
						System.out.println(pathString);
					}
					}
					else{
//						TreePath path = new TreePath(dayNode.getPath());
//						TrajectoryList.add(path);
//						System.out.println(path);

					}

				}
			} else {
//				TreePath path = new TreePath(grNode.getPath());
//				TrajectoryList.add(path);

//				System.out.println(path);

			}
		}
		System.out.println(TrajectoryList);
		return TrajectoryList;
	}

}
