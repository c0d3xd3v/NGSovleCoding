import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

ColumnLayout {
    Layout.bottomMargin: 0
    Layout.margins: 5
    Layout.fillWidth: true


    function updateFunctionList() {
        console.log("update")
        console.log(MainCtrl.getFunctionNames())
        comboBox.model = MainCtrl.getFunctionNames()
    }
    
    Component.onCompleted: function() {
        close()
        MainCtrl.meshLoaded.connect(updateFunctionList)
    }

    RowLayout {
        id: rowLayout1

        //uniformCellSizes: true
        Layout.margins: 0
        Layout.fillHeight: false
        Layout.fillWidth: true

        Label {
            id: label
            opacity: 0.789
            text: qsTr("Gridfunction")
            horizontalAlignment: Text.AlignHCenter
            Layout.fillWidth: false
            Layout.fillHeight: true
            Layout.bottomMargin: 0
            //Layout.fillWidth: true
            Layout.preferredWidth: rowLayout1.width/2
        }

        ComboBox {
            id: comboBox
            opacity: 0.789
            Layout.fillWidth: true
            //Layout.fillWidth: true
            Layout.preferredWidth: rowLayout1.width/2
            onCurrentTextChanged: function() {
                MainCtrl.selectFunctionByName(currentText)
            }
        }

    }



    RowLayout {
        id: rowLayout
        Layout.margins: 0
        Layout.fillHeight: false
        Layout.fillWidth: true

        Label {
            text: qsTr("displacement scale")
            horizontalAlignment: Text.AlignHCenter
            Layout.fillWidth: true
            Layout.preferredWidth: rowLayout.width/2
        }

        TextField {
            id: textField
            Layout.fillWidth: true
            Layout.preferredWidth: rowLayout.width/2
            placeholderText: qsTr("Text Field")
        }
    }


    Item {
        id: item1
        width: 200
        height: 200
        Layout.fillHeight: true
        Layout.preferredWidth: rowLayout1.width

    }
}
